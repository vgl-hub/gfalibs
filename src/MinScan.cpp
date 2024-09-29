/*******************************************************************************************
 *
 *  Package for:
 *    Scan a read buffer of data distributing super-mer streams to buckets.
 *    Read a super-mer stream extracting k-mers for addition to hash tables.
 *
 *  Author:  Gene Myers
 *  Date  :  Sept. 2024
 *
 *******************************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/uio.h>
#include <unistd.h>
#include <fcntl.h>
#include <math.h>
#include <pthread.h>
#include <sys/resource.h>

#include "MinScan.h"

#undef   DEBUG_SETUP
#undef   DEBUG_SCAN
#undef   DEBUG_TRANSMIT
#undef   DEBUG_RECIEVE
#undef   DEBUG_COMPRESSION

typedef unsigned char      uint8;
typedef unsigned long long uint64;

#define BUFFER_LEN  262144  // Size of individual output buffers in uint64's (128 of these)
                            //   => 256MB altogether


/*******************************************************************************************
 *
 *  Initialization (non-threaded) of the package
 *
 *******************************************************************************************/

static char Tran[128] =
  { 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4,
    4, 0, 4, 1, 4, 4, 4, 2,
    4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 3, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4,
    4, 0, 4, 1, 4, 4, 4, 2,
    4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 3, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4,
  };

static uint64  Comp[256];   //  DNA complement of packed byte

static int     Kmer;
static int     Mmer;

static uint64  CHigh[4];
static uint64  NMask;
static int     MaxSuper;
static int     MSWords;
static int     ModMask;
static int     KBits;
static uint64  KMask;
static int     KWords;
static int     SBits;
static uint64  SMask;
static int     PShift;

static pthread_mutex_t BMutex[128];

void Init_Genes_Package(int kmer, int mmer)
{ int i, x;

  Kmer = kmer;
  Mmer = mmer;

  CHigh[0] = 0x3llu<<(2*(mmer-1));
  CHigh[1] = 0x2llu<<(2*(mmer-1));
  CHigh[2] = 0x1llu<<(2*(mmer-1));
  CHigh[3] = 0;
  if (mmer == 32)
    NMask = 0xffffffffffffffffllu;
  else
    NMask = (0x1ll<<(2*mmer))-1;

  PShift = (mmer-7)/2;

  ModMask = 1;
  while (ModMask < kmer)
    ModMask <<= 1;
  ModMask -= 1;

  MaxSuper = kmer-mmer;
  MSWords  = ((Kmer+MaxSuper)-1)/32 + 1;

  KBits  = 2*Kmer;
  KWords = (Kmer-1)/32 + 1;
  KMask = 0xffffffffffffffffllu;
  if ((KBits & 0x3f) > 0)
    KMask -= (0x1llu << (64 - (KBits&0x3f))) - 1;

  SBits = 0;
  x = 1;
  while (x <= MaxSuper)
    { x <<= 1;
      SBits += 1;
    }
  SMask = (0x1llu << SBits)-1;

  for (i = 0; i < 128; i++)
    pthread_mutex_init(BMutex+i,NULL);

  { int l0, l1, l2, l3;   //  Compute byte complement table

    i = 0;
    for (l0 = 3; l0 >= 0; l0 -= 1)
     for (l1 = 12; l1 >= 0; l1 -= 4)
      for (l2 = 48; l2 >= 0; l2 -= 16)
       for (l3 = 192; l3 >= 0; l3 -= 64)
         Comp[i++] = (l3 | l2 | l1 | l0);
  }

#ifdef DEBUG_SETUP
  printf("Kmer     = %d\n",Kmer);
  printf("KBits    = %d\n",KBits);
  printf("KWords   = %d\n",KWords);
  printf("KMask    = %016llx\n",KMask);
  printf("Mmer     = %d\n",Mmer);
  printf("MaxSuper = %d\n",MaxSuper);
  printf("MSWords  = %d\n",MSWords);
  printf("SBits    = %d\n",SBits);
  printf("SMask    = %016llx\n",SMask);
  printf("NMask    = %016llx\n",NMask);
  printf("CHigh[0] = %016llx\n",CHigh[0]);
  printf("CHigh[1] = %016llx\n",CHigh[1]);
  printf("CHigh[2] = %016llx\n",CHigh[2]);
  printf("CHigh[3] = %016llx\n",CHigh[3]);
  printf("ModMask  = %08x\n",ModMask);
  printf("PShift   = %d\n",PShift);
#endif
}


/*******************************************************************************************
 *
 *  Routines, fully re-entrant, for partitioning data into supermers and distribution to
 *  one of 128 files:
 *     Distribution_Bundle *Begin_Distribution(int *fids)
 *     void Distribute_Sequence(char *reads, int rlen, Distribution_Bundle *_bundle)
 *     void End_Distribution(Distribution_Bundle *_bundle)
 *
 *******************************************************************************************/

typedef struct
  { int     fid;  //  file descriptor
    uint64 *ptr;  //  current word ptr
    int     rem;  //  current bit offset in current word
    uint64 *end;  //  end of buffer ptr (< end => room for a supermer)
    uint64 *buf;  //  buffer
  } Packet;

typedef struct
  { Packet  packs[128];  //  distribution buffers
    uint64 *mzr;         //  minimizer queue of size ModMask+1
  } D_Bundle;

  //  Allocate distribution buffers and minimizer queue.
  //  Initialize the buffers.
    
Distribution_Bundle *Begin_Distribution(int *fids)
{ D_Bundle *bundle;
  Packet   *packs;
  int       maxentry;
  int       i;

  bundle = (D_Bundle *) malloc(sizeof(D_Bundle));
  if (bundle == NULL)
    return (NULL);

  packs = bundle->packs;
  packs[0].buf = (uint64 *) malloc(128*BUFFER_LEN*sizeof(uint64));
  if (packs[0].buf == NULL)
    { free(bundle);
      return (NULL);
    }
  maxentry = ((Kmer+MaxSuper)*2+SBits-1)/64 + 1;
  for (i = 0; i < 128; i++)
    { packs[i].buf = packs[0].buf+i*BUFFER_LEN;
      packs[i].ptr = packs[i].buf+1;
      *packs[i].ptr = 0;
      packs[i].rem = 64;
      packs[i].end = packs[i].buf + (BUFFER_LEN-maxentry);
      packs[i].fid = fids[i];
    }

  bundle->mzr = (uint64 *) malloc((ModMask+1)*sizeof(uint64));
  if (bundle->mzr == NULL)
    { free(packs[0].buf);
      free(bundle);
      return (NULL);
    }

#ifdef DEBUG_SETUP
  printf("maxentry = %d\n",maxentry);
#endif

  return ((Distribution_Bundle *) bundle);
}

  //  Flush the distribution buffers and free all working storage.

void End_Distribution(Distribution_Bundle *_bundle)
{ D_Bundle *bundle = (D_Bundle *) _bundle;
  Packet   *packs;
  int       i;

  packs = bundle->packs;
  for (i = 0; i < 128; i++)
    if (packs[i].ptr > packs[i].buf+1 || packs[i].rem < 64)
      { packs[i].buf[0] = ((packs[i].ptr-packs[i].buf)<<6)+(64-packs[i].rem);

#ifdef DEBUG_TRANSMIT
        if (i == 10)
          { printf("\n   Writing Packet of %lld bits",packs[i].buf[0]);
            if (packs[i].rem == 64)
              printf(" in %ld uint64 words\n\n",packs[i].ptr-packs[i].buf);
            else
              printf(" in %ld uint64 words\n\n",(packs[i].ptr-packs[i].buf)+1);
          }
#endif
        pthread_mutex_lock(BMutex+i);
          if (packs[i].rem == 64)
            write(packs[i].fid,packs[i].buf,sizeof(uint64)*(packs[i].ptr-packs[i].buf));
          else
            write(packs[i].fid,packs[i].buf,sizeof(uint64)*((packs[i].ptr-packs[i].buf)+1));
        pthread_mutex_unlock(BMutex+i);
      }

  free(bundle->mzr);
  free(bundle->packs[0].buf);
  free(bundle);
}

//  Stuff overhang length in SBits and advance ptr,pos

static uint64 *stuff_int(uint64 ohang, int *pos, uint64 *ptr)
{ int rem;

  rem = *pos;
  if (rem > SBits)
    { rem -= SBits;
      *ptr |= (ohang << rem);
    }
  else if (rem == SBits)
    { *ptr++ |= ohang;
      rem = 64;
      *ptr = 0;
    }
  else
    { *ptr++ |= (ohang >> (SBits-rem));
      rem += 64 - SBits;
      *ptr = (ohang << rem);
    }
  *pos = rem;
  return (ptr);
}

  //  Stuff seq in 2len bits and advance ptr,pos

static uint64 *stuff_seq(uint8 *seq, int len, int *pos, uint64 *ptr)
{ int    i, rem;
  uint64 v;

  rem = *pos;
  for (i = 0; i < len; i++)
    { v = (uint64) seq[i];
      if (rem > 2)
        { rem -= 2;
          *ptr |= (v << rem);
        }
      else if (rem == 2)
        { *ptr++ |= v;
          rem  = 64;
          *ptr = 0;
        }
      else
        { *ptr++ |= (v >> 1);
          rem  = 64;
          *ptr = (v << rem);
        }
    }
  *pos = rem;
  return (ptr);
}

  //  Transmit supermer seq with overhang ohang to current packet for buffer buck.
  //    Flush the buffer if it becomes full (ptr >= end)

static void transmit(int buck, Packet *pack, uint8 *seq, int ohang)
{ uint64 *ptr;
  int     rem;
#ifdef DEBUG_COMPRESSION
  uint64 *dptr;
#endif

  pack += buck;
  ptr   = pack->ptr;
  rem   = pack->rem;
  
#ifdef DEBUG_TRANSMIT
  if (buck == 10)
    { printf("   Transmit %d %d to bucket %d\n       ",ohang+Kmer,ohang,buck);
      for (int i = 0; i < ohang+Kmer; i++)
        printf("%1d",seq[i]);
      printf("\n");
    }
#endif
#ifdef DEBUG_COMPRESSION
  printf("      %ld:%d\n",ptr-pack->buf,rem);
  dptr = ptr;
#endif

  ptr = stuff_int((uint64) ohang, &rem, ptr);
  ptr = stuff_seq(seq, ohang+Kmer, &rem, ptr);

#ifdef DEBUG_COMPRESSION
  printf("      ");
  while (dptr <= ptr)
    printf(" %016llx",*dptr++);
  printf("\n      %ld:%d\n",ptr-pack->buf,rem);
#endif

  if (ptr >= pack->end)
    { pack->buf[0] = ((ptr-pack->buf)<<6)+(64-rem);

#ifdef DEBUG_TRANSMIT
      if (buck == 10)
        { printf("\n   Writing Packet of %lld bits",pack->buf[0]);
          if (rem == 64)
            printf(" in %ld uint64 words\n\n",ptr-pack->buf);
          else
            printf(" in %ld uint64 words\n\n",(ptr-pack->buf)+1);
        }
#endif
      pthread_mutex_lock(BMutex+buck);
        if (rem == 64)
          write(pack->fid,pack->buf,sizeof(uint64)*(ptr-pack->buf));
        else
          write(pack->fid,pack->buf,sizeof(uint64)*((ptr-pack->buf)+1));
      pthread_mutex_unlock(BMutex+buck);

      ptr = pack->buf+1;
      *ptr = 0;
      rem = 64;
    }

  pack->ptr = ptr;
  pack->rem = rem;
}

  //  Decompose sequences in buffer reads of length rlen into supermers and distribute
  //    to one of 128 files.  The sequences are assumed to be over [acgtACGT] separated
  //    by a single [nN].  Note carefully the last sequences must be followed also by
  //    an [Nn].  The routine does modify the contents of reads mapping, [acgt] to [0123]
  //    to make computing minimizers easier.

void Distribute_Sequence(char *reads, int rlen, Distribution_Bundle *_bundle)
{ D_Bundle *bundle = (D_Bundle *) _bundle;
  Packet *packs    = bundle->packs;
  uint64 *mzr      = bundle->mzr;
  int     kmer1    = Kmer-1;

  int     beg, end, len;
  uint8  *seq, *s, *r;
  int     i, j, last, ohang;
  uint64  x, n, c;
  uint64  mz, mb;
  int     mi;

  seq = (uint8 *) reads;
  if (Tran[seq[rlen-1]] != 4)
    { fprintf(stderr,"String does not end with an N\n");
      exit (1);
    }
  beg = 0;
  for (end = 0; end < rlen; end++)
    { x = seq[end] = Tran[seq[end]];
      if (x < 4)
        continue;
      len = end-beg;
      if (len < Kmer)
        { beg = end+1;
          continue;
        }

      s = seq+beg;
      r = s-kmer1;

      last = kmer1;
      mb = NMask;
      c = n = 0;
      for (i = 0; i < Mmer-1; i++)
        { x = s[i];
          n = (n<<2) | x;
          c = (c>>2) | CHigh[x];
#ifdef DEBUG_SCAN
          printf(" %4d: %1lld %0*llx %0*llx\n",i,x,2*Mmer,n,2*Mmer,c);
#endif
        }
      for (i = Mmer-1; i < Kmer; i++)
        { x = s[i];
          n = ((n<<2) & NMask) | x;
          c = (c>>2) | CHigh[x];
          if (n < c)
            mz = n;
          else
            mz = c;
          mzr[i] = mz;
          if (mz <= mb)
            { mi = i;
              mb = mz;
            }
#ifdef DEBUG_SCAN
          printf(" %4d: %1lld %0*llx %0*llx :: %4d %0*llx\n",i,x,2*Mmer,n,2*Mmer,c,mi,2*Mmer,mb);
#endif
        }
      for (i = Kmer; i < len; i++)
        { x = s[i];
          n = ((n<<2) & NMask) | x;
          c = (c>>2) | CHigh[x];
          if (n < c)
            mz = n;
          else
            mz = c;
          mzr[i&ModMask] = mz;
#ifdef DEBUG_SCAN
          printf(" %4d: %1lld %0*llx %0*llx :: %4d %0*llx\n",i,x,2*Mmer,n,2*Mmer,c,mi,2*Mmer,mb);
#endif
          if (i > mi+MaxSuper)
            { ohang = (i-last)-1;
              transmit((int) (mb&0x7fllu),packs,r+last,ohang); // Super-mer to i-1;
              last = i;
              mb = mzr[(++mi)&ModMask];
              for (j = mi+1; j <= i; j++)
                { mz = mzr[j&ModMask];
                  if (mz <= mb)
                    { mi = j;
                      mb = mz;
                    }
                }
#ifdef DEBUG_SCAN
              printf("%*s:: %4d %0*llx  Forced\n",4*Mmer+11,"",mi,2*Mmer,mb);
#endif
            }
          else if (mz < mb)
            { ohang = (i-last)-1;
              transmit((int) (mb&0x7fllu),packs,r+last,ohang); // Super-mer to i-1;
              last = i;
              mi = i;
              mb = mz;
#ifdef DEBUG_SCAN
              printf("%*s:: %4d %0*llx  New Min\n",4*Mmer+11,"",mi,2*Mmer,mb);
#endif
            }
        }
      ohang = (len-last)-1;
      transmit((int) (mb&0x7fllu),packs,r+last,ohang); // Super-mer to i-1;
  
      beg = end+1;
    }
}


/*******************************************************************************************
 *
 *  Basic routines, fully re-entrant, for processing a bit-packed supermer array:
 *       Scan_Bundle *Begin_Supermer_Scan(uint64 *data, uint64 size)
 *       uint64 *New_Supermer_Buffer()
 *       int Next_Supermer(uint64 *super, Scan_Bundle *_bundle)
 *       void End_Supermer_Scan(Scan_Bundle *_bundle)
 *
 *******************************************************************************************/

typedef struct
  { uint64 *ptr;   // current ptr/rem position in array
    int     rem;
    uint64 *eptr;  // last position in the current packet
    int     erem;
    uint64 *dend;  // ptr to the end of the data array
  } S_Bundle;

  //  Setup to scan/process a section (or all) of a data array starting
  //    at data of length size uint64 words.

static uint64 *DATA;

Scan_Bundle *Begin_Supermer_Scan(uint64 *data, uint64 size)
{ S_Bundle *bundle;
  uint64    totbits;

  bundle = (S_Bundle *) malloc(sizeof(S_Bundle));
  if (bundle == NULL)
    return (NULL);

  totbits = *data;
  bundle->erem = 64 - (totbits & 0x3f);
  bundle->eptr = data+(totbits>>6);
  bundle->rem  = 64;
  bundle->ptr  = data+1;
  bundle->dend = data+size;

DATA = bundle->ptr;

#ifdef DEBUG_RECIEVE
  printf("\n   Scanning Packet of %lld bits",totbits);
  if (bundle->erem == 64)
    printf(" covering %lld uint64 words\n\n",totbits>>6);
  else
    printf(" covering %lld uint64 words\n\n",(totbits>>6)+1);
  fflush(stdout);
#endif

  return ((Scan_Bundle *) bundle);
}

  // End a scan by free the bundle

void End_Supermer_Scan(Scan_Bundle *_bundle)
{ S_Bundle *bundle = (S_Bundle *) _bundle;

  free(bundle);
  return;
}

  // Allocate a return an array big enough to hold any supermer
  //   Return NULL if there was not enough memory (unlikely).
  //   The user is responsible for ultimately free'ing the returned array.

uint64 *New_Supermer_Buffer()
{ uint64 *super;

  super = (uint64 *) malloc(sizeof(uint64)*(Kmer+MaxSuper));
  if (super == NULL)
    return (NULL);
  return (super);
}

  // Assumes ptr/pos is just before an overhang length.  Unpack the next SBits into
  //   ohang and advance ptr/pos over it returning ptr.

static uint64 *unstuff_int(int *ohang, int *pos, uint64 *ptr)
{ int rem;

  rem = *pos;
  if (SBits < rem)
    { rem -= SBits;
      *ohang = (((*ptr) >> rem) & SMask);
    }
  else if (SBits == rem)
    { *ohang = (*ptr++) & SMask;
      rem  = 64;
    }
  else
    { uint64 x = (*ptr++) << (SBits-rem);
      rem += 64 - SBits;
      *ohang = (x | (*ptr >> rem)) & SMask;
    }
  *pos = rem;
  return (ptr);
}

  // Assumes ptr/pos is just before a 2-bit packet supermer of length leng.
  //   Unpack it into the array super and advance ptr/pos over it return ptr.

static uint64 *unstuff_seq(uint64 *s, int len, int *pos, uint64 *ptr)
{ int    i, rem;
  uint64 x;

  rem = *pos;
  for (i = 32; i <= len; i += 32)
    if (rem == 64)
      *s++ = *ptr++;
    else
      { x = ((*ptr++) << (64-rem));
        *s++ = (x | (*ptr >> rem));
      }
  i = (len & 0x1f) << 1;
  if (i > 0)
    { if (rem > i)
        { rem -= i;
          *s++ = ((*ptr >> rem) << (64-i));
        }
      else if (rem == i)
        { *s++ = *ptr++ << (64-rem);
          rem = 64;
        }
      else
        { x = ((*ptr++) << (64-rem));
          rem += 64-i;
          *s++ = (x | (*ptr >> rem) << (64-i));
        }
    }
  *pos = rem;
  return (ptr);
}

  //  Assumes ptr/pos is just before an overhang length follwed by a 2-bit
  //    packed supermer of the appropriate length. Get the next supermer
  //    packed in super and return its length.  0 is returned if the end
  //    of the data section is reached.  The bundle pointer is advanced
  //    over the overhang and supermer bits, processing packet boundaries
  //    as necessary.

int Next_Supermer(uint64 *super, Scan_Bundle *_bundle)
{ S_Bundle *bundle = (S_Bundle *) _bundle;
  uint64 totbits;
  int    ohang, len;

  if (bundle->ptr >= bundle->eptr && bundle->rem <= bundle->erem)
    { if (bundle->rem < 64)
        bundle->ptr += 1;
      if (bundle->ptr >= bundle->dend)
        return (0);
      totbits = *bundle->ptr;
      bundle->erem = 64 - (totbits & 0x3f);
      bundle->eptr = bundle->ptr+(totbits>>6);
      bundle->rem  = 64;
      bundle->ptr += 1;
#ifdef DEBUG_RECIEVE
      printf("\n   Scanning Packet of %lld bits",totbits);
      if (bundle->erem == 64)
        printf(" covering %lld uint64 words\n\n",totbits>>6);
      else
        printf(" covering %lld uint64 words\n\n",(totbits>>6)+1);
      fflush(stdout);
#endif
    }

  bundle->ptr = unstuff_int(&ohang,&bundle->rem,bundle->ptr);
  len = ohang+Kmer;
  bundle->ptr = unstuff_seq(super,len,&bundle->rem,bundle->ptr);

#ifdef DEBUG_RECIEVE
  { int i, w, x;

    printf("   Transmit %d %d to bucket 10\n       ",len,ohang);
    x = 62;
    w = 0;
    for (i = 0; i < len; i++)
      { printf("%1lld",(super[w]>>x)&0x3llu);
        if (x == 0)
          { w += 1;
            x = 64;
          }
        x -= 2;
      }
    printf("\n");
    fflush(stdout);
  }
#endif

  return (len);
}


/*******************************************************************************************
 *
 *  Routines, fully re-entrant, for traversing a bit-packed supermer array either packet
 *  by packet, or super-mer by super-mer:
 *     uint64 *Skip_To_Next_Packet(Scan_Bundle *_bundle)
 *     int Get_Kmer_Count(Scan_Bundle *_bundle)
 *     void Skip_Kmers(int len, Scan_Bundle *_bundle)
 *
 *******************************************************************************************/

  //  Given an intialized scan bundle at any point in a data segment, advance to the start of
  //    the next packet (always on a uint64 word boundary) and return a pointer to the start
  //    of said packet.  If there is no next packet then return NULL.

uint64 *Skip_To_Next_Packet(Scan_Bundle *_bundle)
{ S_Bundle *bundle = (S_Bundle *) _bundle;
  uint64    totbits;

  bundle->ptr = bundle->eptr;
  if (bundle->erem < 64)
    bundle->ptr += 1;
  if (bundle->ptr >= bundle->dend)
    return (NULL);
  totbits = *bundle->ptr;
  bundle->erem = 64 - (totbits & 0x3f);
  bundle->eptr = bundle->ptr+(totbits>>6);
  bundle->rem  = 64;
  bundle->ptr += 1;
  return (bundle->ptr-1);
}

  //  Assumes the bundle points at the start of the overhang length of a supermer record.
  //  Get and scan past the overhang length and return the # of k-mers in the supermer,
  //  i.e. over hang length + 1;

int Get_Kmer_Count(Scan_Bundle *_bundle)
{ S_Bundle *bundle = (S_Bundle *) _bundle;
  uint64 totbits;
  int    ohang;

  if (bundle->ptr >= bundle->eptr && bundle->rem <= bundle->erem)
    { if (bundle->rem < 64)
        bundle->ptr += 1;
      if (bundle->ptr >= bundle->dend)
        return (0);
      totbits = *bundle->ptr;
      bundle->erem = 64 - (totbits & 0x3f);
      bundle->eptr = bundle->ptr+(totbits>>6);
      bundle->rem  = 64;
      bundle->ptr += 1;
#ifdef DEBUG_RECIEVE
      printf("\n   Scanning Packet of %lld bits",totbits);
      if (bundle->erem == 64)
        printf(" covering %lld uint64 words\n\n",totbits>>6);
      else
        printf(" covering %lld uint64 words\n\n",(totbits>>6)+1);
      fflush(stdout);
#endif
    }

  bundle->ptr = unstuff_int(&ohang,&bundle->rem,bundle->ptr);

#ifdef DEBUG_RECIEVE
  printf("   Transmit %d %d to bucket 10\n       ",ohang+Kmer,ohang);
  fflush(stdout);
#endif

  return (ohang+1);
}

  //  Assumes the bundle points at the start of a 2-bit supermer containing count k-mer.
  //  Advance the bundle position over the supermer to the start of the next supermer record.

void Skip_Kmers(int count, Scan_Bundle *_bundle)
{ S_Bundle *bundle = (S_Bundle *) _bundle;
  int rem;

  count += Kmer-1;
  bundle->ptr += (count>>5);
  count = (count & 0x1f) << 1;
  rem   = bundle->rem;
  if (rem > count)
    rem -= count;
  else
    { rem += 64-count;
      bundle->ptr += 1;
    }
  bundle->rem = rem;
}


/*******************************************************************************************
 *
 *  Routines, fully re-entrant, for fetching hashes and k-mers from a packed supermer array
 *  given an 'offset' into the array:
 *     uint64 Current_Offset(uint64 *finger, Scan_Bundle *_bundle)
 *     void Get_Hash(int *dir, uint64 *hash, uint64 *finger, uint64 offset)
 *     void Get_Canonical_Kmer(uint64 *super, int *dir, uint64 hash, uint64 *finger, uint64 offset)
 *
 *******************************************************************************************/

  //  Return a uint64 that precisely encodes the position of the current bundle with respect to
  //    a "finger" pointer into the data array.  The bundle must be at a position beyond the
  //    finger so that the "offset" is positive.  Arithmetic with offsets works, e.g. offset+2
  //    is 2 bits further from the finger.

uint64 Current_Offset(uint64 *finger, Scan_Bundle *_bundle)
{ S_Bundle *bundle = (S_Bundle *) _bundle;
  uint64    offset;

  offset = bundle->ptr-finger;
  offset = (offset << 6) | (64-bundle->rem);
  return (offset);
}

  // Get the complement of the 32bp *ending* at position ptr/rem

static inline uint64 get_comp64(int rem, uint64 *ptr)
{ uint64 x, w;
  int    i;

  if (rem == 64)
    w = ptr[-1];
  else
    w = (*ptr >> rem) | (ptr[-1] << (64-rem));
  x = (Comp[w&0xffllu]<<56);
  for (i = 8; i < 64; i += 8)
    x |= (Comp[(w >> i) & 0xffllu] << (56-i));
  return (x);
}

  // Get the 32bp beginning at position ptr/rem

static inline uint64 get_norm64(int rem, uint64 *ptr)
{ uint64 x;

  if (rem == 64)
    x = *ptr;
  else
    x = (*ptr << (64-rem)) | (ptr[1] >> rem);
  return (x);
}

  //  Get the hash (1st 32bp or k-mer + zero-padding if k < 32) of the canonical k-mer at the
  //    given offset from finger.  Return +1 if the hash is from the forward direction,
  //    -1 if from the complement direction, and 0 if the two are equal.

int Get_Hash(uint64 *hash, uint64 *finger, uint64 offset)
{ uint64 *ptr;
  int     rem;
  uint64  n, c;

  rem = 64-(offset & 0x3f);
  ptr = finger + (offset>>6);
  n = get_norm64(rem,ptr);
  if (rem < KBits)
    { rem  = KBits-rem;
      ptr += (rem>>6) + 1;
      rem  = 64 - (rem & 0x3f);
    }
  else
    rem -= KBits;
  c = get_comp64(rem,ptr);

  if (Kmer < 32)
    { c &= KMask;
      n &= KMask;
    }
  if (c < n)
    { *hash = c;
      return (-1);
    }
  else if (c == n)
    { *hash = c;
      return (0);
    }
  else
    { *hash = n;
      return (1);
    }
}

  //  Get the canonical k-mer at the offset from finger whose hash is hash that comes from the end
  //    specified by dir (both obtained from Get_Hash).  1 is returned if the k-mer is in the
  //    forward direction and -1 if in the complement direction.

int Get_Canonical_Kmer(uint64 *super, int dir, uint64 hash, uint64 *finger, uint64 offset)
{ uint64 *ptr, *ctr;
  uint64  n, c;
  int     rem, cem;
  int     i;

  *super++ = hash;
  if (Kmer <= 32)
    { if (dir >= 0)
        return (1);
      else
        return (-1);
    }

  rem = (offset & 0x3f)+1;
  ptr = finger + (offset>>6);
  if (dir <= 0)
    { cem  = KBits-rem;
      ctr  = ptr + (cem>>6);
      cem  = 64 - (cem & 0x3f);
    }

  i = 0;
  if (dir == 0)
    { for (i = 0; i < KWords; i++)
        { n = get_norm64(rem,ptr++);
          c = get_comp64(cem,ctr--);
          if (c < n)
            { dir = -1;
              *super++ = c;
              break;
            }
          if (c > n)
            { dir = 1;
              *super++ = n;
              break;
            }
          *super++ = c;
         }
       if (i >= KWords)
         { n = get_norm64(rem,ptr) & KMask;
           c = get_comp64(cem,ctr) & KMask;
           if (c < n)
             { *super = c;
               return (-1);
             }
           else
             { *super = n;
               return (1);
             }
         }
       i += 1;
    }

  if (dir < 0)
    { for ( ; i < KWords; i++)
        *super++ = get_comp64(cem,ctr--);
      *super = get_comp64(cem,ctr) & KMask;
      return (-1);
    }

  for ( ; i < KWords; i++)
    *super++ = get_norm64(rem,ptr++);
  *super = get_norm64(rem,ptr) & KMask;
  return (1);
}
