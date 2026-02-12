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

#if defined(_WIN32) || defined(_WIN64)
// MinGW/Windows doesn't ship <sys/uio.h>.
// If you only need iovec, define it here.
#include <cstddef>   // size_t
struct iovec {
	void*  iov_base;
	size_t iov_len;
};

// If MinScan actually calls readv/writev, add a small wrapper
// (or refactor to plain read()/write()).
#else
#include <sys/uio.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#include <fcntl.h>
#include <math.h>
#include <pthread.h>
#include <sys/resource.h>

#include "MinScan.h"

#undef    DEBUG_SETUP
#undef    DEBUG_SCAN
#undef    DEBUG_TRANSMIT
#undef    DEBUG_RECIEVE
#undef    DEBUG_COMPRESSION

#if defined(DEBUG_SETUP) || defined(DEBUG_SCAN)

static char dna[4] = { 'a', 'c', 'g', 't' };

static char *fmer[256], _fmer[1280];

static void setup_fmer_table()
{ char *t;
  int   i, l3, l2, l1, l0;

  i = 0;
  t = _fmer;
  for (l3 = 0; l3 < 4; l3++)
   for (l2 = 0; l2 < 4; l2++)
	for (l1 = 0; l1 < 4; l1++)
	 for (l0 = 0; l0 < 4; l0++)
	   { fmer[i] = t;
		 *t++ = dna[l3];
		 *t++ = dna[l2];
		 *t++ = dna[l1];
		 *t++ = dna[l0];
		 *t++ = 0;
		 i += 1;
	   }
}

#endif


/*******************************************************************************************
 *
 *  Initialization and training (non-threaded) of the package
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

static char Invert[5] = { 'a', 'c', 'g', 't', 'n' };

static uint64  Comp[256];   //  DNA complement of packed byte

static int     Kmer;
static int     Mmer;

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

void Init_Genes_Package(int kmer, int mmer)
{ int i, x;

  Kmer = kmer;
  Mmer = mmer;

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

  { int l0, l1, l2, l3;   //  Compute byte complement table

	i = 0;
	for (l0 = 3; l0 >= 0; l0 -= 1)
	 for (l1 = 12; l1 >= 0; l1 -= 4)
	  for (l2 = 48; l2 >= 0; l2 -= 16)
	   for (l3 = 192; l3 >= 0; l3 -= 64)
		 Comp[i++] = (l3 | l2 | l1 | l0);
  }

#ifdef DEBUG_SETUP
  printf("\n");
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
  printf("ModMask  = %08x\n",ModMask);
  printf("PShift   = %d\n",PShift);
#endif
}


/*******************************************************************************************
 *
 *  Training phase: a single call to Train_Genes_Package
 *
 *******************************************************************************************/

static int              Fnum;        // # of buckets
static pthread_mutex_t *BMutex;      // Array of Fnum Mutex's for each bucket's IO
static int              Buffer_Len;  // Size of individual output buffers in uint64's
									 //   so that 1GB is used altogether.
static uint64           TMap[256];   // Map 4-tuples to code/value
static uint64           CMap[256];   // Map 4-tuple code to code of complement high-order part
static int64           *Count;       // Trimmed bucket mapping tree

static int FREQ[256];

static int FSORT(const void *l, const void *r)
{ int x = *((int *) l);
  int y = *((int *) r);
  return (FREQ[x] - FREQ[y]);
}

#undef DEBUG_SCAN1

static int Test_Mmer_Choice(char *reads, int rlen, int64 *all)
{ int     kmer1 = Kmer-1;

  int     beg, end, len;
  uint8  *seq, *s;
  int     i, j, last, ohang;
  uint64  x, n, c, N[4], C[4];
  uint64  mz, mb, *mzr;
  int     mi;
  int64   hzero, hall;

  mzr = (uint64 *) malloc((ModMask+1)*sizeof(uint64));

#ifdef DEBUG_SCAN1
  setup_fmer_table();
#endif

  seq = (uint8 *) reads;
  hzero = hall = 0;
  beg = 0;
  for (end = 0; end < rlen; end++)
	{ x = seq[end];
	  if (x < 4)
		continue;
	  len = end-beg;
	  if (len < Kmer)
		{ beg = end+1;
		  continue;
		}

	  s = seq+beg;

	  last = kmer1;
	  mb = mi = NMask;
	  for (j = 0; j < 4; j++)
		N[j] = C[j] = 0;
	  x = (s[0] << 4) | (s[1] << 2) | s[2];
	  for (i = 3; i < Mmer-1; i++)
		{ x = ((x<<2) | s[i]) & 0xff;
		  j = (i & 0x3);
		  N[j] = (N[j]<<8) | TMap[x];
		  C[j] = (C[j]>>8) | CMap[x];
#ifdef DEBUG_SCAN1
		  printf(" %4d: %s %0*llx %0*llx\n",i,fmer[x],Mmer/2,N[j],Mmer/2,C[j]);
#endif
		}
	  for (i = Mmer-1; i < Kmer; i++)
		{ x = ((x<<2) | s[i]) & 0xff;
		  j = (i & 0x3);
		  n = N[j] = ((N[j]<<8) & NMask) | TMap[x];
		  c = C[j] = (C[j]>>8) | CMap[x];
#ifdef DEBUG_SCAN1
		  printf(" %4d: %s %0*llx %0*llx\n",i,fmer[x],Mmer/2,n,Mmer/2,c);
#endif
		  if (n < c)
			mz = n;
		  else
			mz = c;
		  mzr[i] = mz;
		  if (mz <= mb)
			{ mi = i;
			  mb = mz;
			}
		}
	  for (i = Kmer; i < len; i++)
		{ x = ((x<<2) | s[i]) & 0xff;
		  j = (i & 0x3);
		  n = N[j] = ((N[j]<<8) & NMask) | TMap[x];
		  c = C[j] = (C[j]>>8) | CMap[x];
#ifdef DEBUG_SCAN1
		  printf(" %4d: %s %0*llx %0*llx :: %4d %0*llx\n",i,fmer[x],Mmer/2,n,Mmer/2,c,mi,Mmer/2,mb);
#endif
		  if (n < c)
			mz = n;
		  else
			mz = c;
		  mzr[i&ModMask] = mz;
#ifdef DEBUG_SCAN1
		  printf(" %4d: %s %0*llx %0*llx :: %4d %0*llx\n",i,fmer[x],Mmer/2,n,Mmer/2,c,mi,Mmer/2,mb);
#endif
		  if (i > mi+MaxSuper)
			{ ohang = (i-last)-1;
			  if (mb == 0)
				hzero += Kmer+ohang;
			  hall += Kmer+ohang;
			  last = i;

			  mb = NMask;
			  j  = mi+1;
			  for (mi = -1; j <= i; j++)
				{ mz = mzr[j&ModMask];
				  if (mz <= mb)
					{ mi = j;
					  mb = mz;
					}
				}
#ifdef DEBUG_SCAN1
			  printf("%*s:: %4d %0*llx  Forced\n",4*Mmer+11,"",mi,Mmer/2,mb);
#endif
			}
		  else if (mz < mb)
			{ ohang = (i-last)-1;
			  if (mb == 0)
				hzero += Kmer+ohang;
			  hall += Kmer+ohang;
			  last = i;
			  mi = i;
			  mb = mz;
#ifdef DEBUG_SCAN1
			  printf("%*s:: %4d %0*llx  New Min\n",4*Mmer+11,"",mi,Mmer/2,mb);
#endif
			}
		}
	  ohang = (len-last)-1;
	  if (mb == 0)
		hzero += Kmer+ohang;
	  hall += Kmer+ohang;
  
	  beg = end+1;
	}

  free(mzr);

#ifdef DEBUG_SETUP
  printf("\nTesting mmer = %d: %lld / %lld = %d\n",Mmer,hall,hzero,(int) (hall/hzero));
#endif

  *all = hall;
  return ((int) (hall/hzero));
}

#undef DEBUG_SCAN2

static void Count_Map_Tree(char *reads, int rlen, int64 *count)
{ int     kmer1 = Kmer-1;
  int     pad   = 2*(Mmer-5);

  int     beg, end, len;
  uint8  *seq, *s;
  int     i, j, last, ohang;
  uint64  x, n, c, N[4], C[4];
  uint64  mz, mb, *mzr;
  int     mi;
  int     y, p, v;

  mzr = (uint64 *) malloc((ModMask+1)*sizeof(uint64));

#ifdef DEBUG_SCAN2
  setup_fmer_table();
#endif

  seq = (uint8 *) reads;
  beg = 0;
  for (end = 0; end < rlen; end++)
	{ x = seq[end];
	  if (x < 4)
		continue;
	  len = end-beg;
	  if (len < Kmer)
		{ beg = end+1;
		  continue;
		}

	  s = seq+beg;

	  last = kmer1;
	  mb = mi = NMask;
	  for (j = 0; j < 4; j++)
		N[j] = C[j] = 0;
	  x = (s[0] << 4) | (s[1] << 2) | s[2];
	  for (i = 3; i < Mmer-1; i++)
		{ x = ((x<<2) | s[i]) & 0xff;
		  j = (i & 0x3);
		  N[j] = (N[j]<<8) | TMap[x];
		  C[j] = (C[j]>>8) | CMap[x];
#ifdef DEBUG_SCAN2
		  printf(" %4d: %s %0*llx %0*llx\n",i,fmer[x],Mmer/2,N[j],Mmer/2,C[j]);
#endif
		}
	  for (i = Mmer-1; i < Kmer; i++)
		{ x = ((x<<2) | s[i]) & 0xff;
		  j = (i & 0x3);
		  n = N[j] = ((N[j]<<8) & NMask) | TMap[x];
		  c = C[j] = (C[j]>>8) | CMap[x];
#ifdef DEBUG_SCAN2
		  printf(" %4d: %s %0*llx %0*llx\n",i,fmer[x],Mmer/2,n,Mmer/2,c);
#endif
		  if (n < c)
			mz = n;
		  else
			mz = c;
		  mzr[i] = mz;
		  if (mz <= mb)
			{ mi = i;
			  mb = mz;
			}
		}
	  for (i = Kmer; i < len; i++)
		{ x = ((x<<2) | s[i]) & 0xff;
		  j = (i & 0x3);
		  n = N[j] = ((N[j]<<8) & NMask) | TMap[x];
		  c = C[j] = (C[j]>>8) | CMap[x];
#ifdef DEBUG_SCAN2
		  printf(" %4d: %s %0*llx %0*llx :: %4d %0*llx\n",i,fmer[x],Mmer/2,n,Mmer/2,c,mi,Mmer/2,mb);
#endif
		  if (n < c)
			mz = n;
		  else
			mz = c;
		  mzr[i&ModMask] = mz;
		  if (i > mi+MaxSuper)
			{ ohang = (i-last)-1;

			  p = pad;
			  y = (mb >> p);
			  v = count[y];
			  while (v < 0)
				{ p -= 2;
				  y  = ((mb >> p) & 0x3) - v;
				  v  = count[y];
				}
			  count[y] += Kmer+ohang;

			  last = i;

			  mb = NMask;
			  j  = mi+1;
			  for (mi = -1; j <= i; j++)
				{ mz = mzr[j&ModMask];
				  if (mz <= mb)
					{ mi = j;
					  mb = mz;
					}
				}
#ifdef DEBUG_SCAN2
			  printf("%*s:: %4d %0*llx  Forced\n",4*Mmer+11,"",mi,Mmer/2,mb);
#endif
			}
		  else if (mz < mb)
			{ ohang = (i-last)-1;

			  p = pad;
			  y = (mb >> p);
			  v = count[y];
			  while (v < 0)
				{ p -= 2;
				  y  = ((mb >> p) & 0x3) - v;
				  v  = count[y];
				}
			  count[y] += Kmer+ohang;

			  last = i;
			  mi = i;
			  mb = mz;
#ifdef DEBUG_SCAN2
			  printf("%*s:: %4d %0*llx  New Min\n",4*Mmer+11,"",mi,Mmer/2,mb);
#endif
			}
		}
	  ohang = (len-last)-1;

	  p = pad;
	  y = (mb >> p);
	  v = count[y];
	  while (v < 0)
		{ p -= 2;
		  y  = ((mb >> p) & 0x3) - v;
		  v  = count[y];
		}
	  count[y] += Kmer+ohang;
  
	  beg = end+1;
	}

  free(mzr);
}

#ifdef DEBUG_SETUP

static void _print_count_tree(int lev, int i, int64 ktot, int64 *count)
{ int j, a;

  if (count[i] >= 0)
	printf(" %10lld %.3f%% (%d)\n",count[i],(100.*count[i])/ktot,lev);
  else
	{ printf("\n");
	  j = -count[i];
	  for (a = 0; a < 4; a++)
		{ printf("%*s -> %d:",2*lev,"",a);
		  _print_count_tree(lev+1,j+a,ktot,count);
		}
	}
}

static void print_count_tree(int64 ktot, int64 *count)
{ int i;

  printf("\nCount Tree:\n");
  for (i = 0; i < 1024; i++)
	{ printf(" %5d:",i);
	  _print_count_tree(0,i,ktot,count);
	}
  fflush(stdout);
}

static void _print_assign_tree(int lev, int i, int64 *count)
{ int j, a;

  if (count[i] >= 0)
	printf(" %10lld (%d)\n",count[i],lev);
  else
	{ printf("\n");
	  j = -count[i];
	  for (a = 0; a < 4; a++)
		{ printf("%*s -> %d:",2*lev,"",a);
		  _print_assign_tree(lev+1,j+a,count);
		}
	}
}

static void print_assign_tree(int64 *count)
{ int i;

  printf("\nAssign Tree:\n");
  for (i = 0; i < 1024; i++)
	{ printf(" %5d:",i);
	  _print_assign_tree(0,i,count);
	}
  fflush(stdout);
}

#endif

static int CSORT(const void *l, const void *r)
{ int64 x = Count[*((int*) l)];
  int64 y = Count[*((int*) r)];
  
  if (x > y)
	return (-1);
  else
	return (1);
}

static void Map_Assignment(int csize, int64 *count, int64 cthresh)
{ int64 *buck;
  int   *perm;
  int    i, j, n, x;
  int64  p;

  perm = (int *) malloc(sizeof(int)*csize);
  buck = (int64 *) malloc(sizeof(int64)*Fnum);

  for (i = 0; i < csize; i++)
	perm[i] = i;

  qsort(perm,csize,sizeof(int),CSORT);

#ifdef DEBUG_SETUP
  printf("\nPieces %d:\n",csize);
  for (i = 0; i < csize; i++)
	if (count[perm[i]] > 0)
	  printf("  %5d: %7lld -> %5.2f%%\n",perm[i],count[perm[i]],(100.*count[perm[i]])/cthresh);
  fflush(stdout);
#endif

  for (i = 0; i < Fnum; i++)
	{ x = perm[i];
	  buck[i]  = count[x];
	  count[x] = i;
	}
  for (i = Fnum; i < csize; i++)
	{ x = perm[i];
	  p = count[x];
	  if (p < 0)      //  only interior nodes remain
		break;
	  if (p == 0)
		{ count[x] = (int) (drand48() * Fnum);   //  never seen?: place at random
		  continue;
		}
	  for (j = 0; j < Fnum; j++)       // place in first bucket it fits in
		if (buck[j] + p <= cthresh)
		  { buck[j] += p;
			count[x] = j;
			break;
		  }
	  if (j >= Fnum)    //  if none, then place in the emptiest bucket
		{ n = 0;
		  for (j = 1; j < Fnum; j++)
			if (buck[j] < buck[n])
			  n = j;
		  buck[n] += p;
		  count[x] = n;
		}
	}

#ifdef DEBUG_SETUP
  { int64 bmax, bmin;

	printf("\nPacking:\n");
	bmax = bmin = buck[0];
	for (i = 0; i < Fnum; i++)
	  { printf("  %7llu -> %5.2f%%\n",buck[i],(100.*buck[i])/cthresh);
		if (bmax < buck[i])
		  bmax = buck[i];
		else if (bmin > buck[i])
		  bmin = buck[i];
	  }
	printf("  Range %lld - %lld => %g%% diff\n",bmin,bmax,(100.*(bmax-bmin))/cthresh);
	fflush(stdout);
  }
#endif

  free(buck);
  free(perm);
}

  //  To start tell the package the kmer size and the number of files the supermers should be
  //    distributed over.  Also give it an initial cache of data of length clen to train
  //    on.  Based on this segment the routine will select a minimizer size to use and develop
  //    a mapping scheme, as per the one used in FastK.  The routine returns 0 if it could
  //    not produce a scheme (should never happen basically) and otherwise the size of the
  //    minimizer it will use.

int Train_Genes_Package(int kmer, int fnum, char *cache, int clen)
{ int     n, x;
  uint8  *seq;
  int     fmr, fln;
  int     perm[256];
  int     csize;
  int64   cthresh, total;

  Fnum       = fnum;
  Buffer_Len = 0x40000000llu / (fnum*sizeof(uint64));
  BMutex     = (pthread_mutex_t *) malloc(sizeof(pthread_mutex_t)*fnum);
  for (n = 0; n < fnum; n++)
	pthread_mutex_init(BMutex+n,NULL);

  for (n = 0; n < 128; n++)
	FREQ[n] = 0;

  seq = (uint8 *) cache;
  fmr = 0;
  fln = 0;
  for (n = 0; n < clen; n++)
	{ x = seq[n] = Tran[seq[n]];
	  if (x == 4)
		{ fmr = 0;
		  fln = 0;
		}
	  else
		{ fmr = ((fmr << 2) | x) & 0xff;
		  fln += 1;
		  if (fln >= 4)
			FREQ[fmr] += 1;
		}
	}

  for (n = 0; n < 256; n++)
	perm[n] = n;

  qsort(perm,256,sizeof(int),FSORT);

#ifdef DEBUG_SETUP
  setup_fmer_table();

  printf("\n4-mer training frequencies in order\n");
  for (n = 0; n < 256; n++)
	printf(" %02x: %s: %d\n",n,fmer[perm[n]],FREQ[perm[n]]);
#endif

  for (n = 0; n < 256; n++)
	TMap[perm[n]] = n;

#ifdef DEBUG_SETUP
  printf("\n4-mer forward map\n");
  for (n = 0; n < 256; n++)
	printf(" %s -> %02llx\n",fmer[n],TMap[n]);
#endif

  Init_Genes_Package(kmer,8);
  for (n = 0; n < 256; n++)
	CMap[n] = (TMap[Comp[n]] << 8);

#ifdef DEBUG_SETUP
  printf("\n4-mer complement map\n");
  for (n = 0; n < 256; n++)
	printf(" %s -> %016llx\n",fmer[n],CMap[n]);
#endif

  for (n = 8; n < 32; n += 4)
	{ if (Test_Mmer_Choice(cache,clen,&total) >= 2*Fnum)
		break;
	  Init_Genes_Package(kmer,n+4);
	  for (x = 0; x < 256; x++)
		CMap[x] = (TMap[Comp[x]] << (2*n));
	}
  if (n >= 32)
	return (0);

  csize   = 1024;
  Count   = (int64 *) malloc(sizeof(int64)*csize);
  cthresh = total / (2*Fnum);

#ifdef DEBUG_SETUP
  printf("Threshold = %lld\n",cthresh);
#endif

  for (n = 0; n < csize; n++)
	Count[n] = 0;

  while (1)
	{ int expand;

	  Count_Map_Tree(cache,clen,Count);

#ifdef DEBUG_SETUP
	  print_count_tree(total,Count);
#endif
 
	  expand = csize;
	  for (n = 0; n < csize; n++)
		if (Count[n] > cthresh)
		  expand += 4;

	  if (expand == csize)
		break;

	  Count = (int64 *) realloc(Count,sizeof(int64)*expand);

	  expand = csize;
	  for (n = 0; n < csize; n++)
		{ if (Count[n] < 0)
			continue;
		  if (Count[n] > cthresh)
			{ Count[n] = -expand;
			  for (x = 0; x < 4; x++)
				Count[expand++] = 0;
			}
		  else
			Count[n] = 0;
		}
	  csize = expand;
	}

  Map_Assignment(csize,Count,total/Fnum);

#ifdef DEBUG_SETUP
  print_assign_tree(Count);
#endif

  for (n = 0; n < clen; n++)
	cache[n] = Invert[seq[n]];

  return (Mmer);
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
  { Packet *packs;       //  distribution buffers
	uint64 *mzr;         //  minimizer queue of size ModMask+1
  } D_Bundle;

  //  Allocate distribution buffers and minimizer queue.
  //    Initialize the buffers.
	
Distribution_Bundle *Begin_Distribution(int *fids)
{ D_Bundle *bundle;
  Packet   *packs;
  int       maxentry;
  int       i;

  bundle = (D_Bundle *) malloc(sizeof(D_Bundle));
  packs  = (Packet *) malloc(sizeof(Packet)*Fnum);
  if (bundle == NULL || packs == NULL)
	return (NULL);
  bundle->packs = packs;

  packs[0].buf = (uint64 *) malloc(Fnum*Buffer_Len*sizeof(uint64));
  if (packs[0].buf == NULL)
	{ free(bundle);
	  return (NULL);
	}
  maxentry = ((Kmer+MaxSuper)*2+SBits-1)/64 + 1;
  for (i = 0; i < Fnum; i++)
	{ packs[i].buf = packs[0].buf+i*Buffer_Len;
	  packs[i].ptr = packs[i].buf+1;
	  *packs[i].ptr = 0;
	  packs[i].rem = 64;
	  packs[i].end = packs[i].buf + (Buffer_Len-maxentry);
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
  for (i = 0; i < Fnum; i++)
	if (packs[i].ptr > packs[i].buf+1 || packs[i].rem < 64)
	  { packs[i].buf[0] = ((packs[i].ptr-packs[i].buf)<<6)+(64-packs[i].rem);

#ifdef DEBUG_TRANSMIT
		if (i == 1)
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
  free(bundle->packs);
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
		  rem  = 63;
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
  if (buck == 1)
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
	  if (buck == 1)
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
  //    to one of Fnum files.  The sequences are assumed to be over [acgtACGT] separated
  //    by a single [nN].  Note carefully the last sequences must be followed also by
  //    an [Nn].  The routine does modify the contents of reads mapping, [acgt] to [0123]
  //    to make computing minimizers easier.

void Distribute_Sequence(char *reads, int rlen, Distribution_Bundle *_bundle)
{ D_Bundle *bundle = (D_Bundle *) _bundle;
  Packet *packs    = bundle->packs;
  uint64 *mzr      = bundle->mzr;
  int     kmer1    = Kmer-1;
  int     pad      = 2*(Mmer-5);

  int     beg, end, len;
  uint8  *seq, *s, *r;
  int     i, j, last, ohang;
  uint64  x, n, c, N[4], C[4];
  uint64  mz, mb;
  int     mi;
  int     y, p, v;

#ifdef DEBUG_SCAN
  setup_fmer_table();
  print_assign_tree(Count);
#endif

  seq = (uint8 *) reads;
  if (Tran[seq[rlen-1]] != 4)
	{ fprintf(stderr,"Read buffer does not end with an N\n");
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
	  mb = mi = NMask;
	  for (j = 0; j < 4; j++)
		N[j] = C[j] = 0;
	  x = (s[0] << 4) | (s[1] << 2) | s[2];
	  for (i = 3; i < Mmer-1; i++)
		{ x = ((x<<2) | s[i]) & 0xff;
		  j = (i & 0x3);
		  N[j] = (N[j]<<8) | TMap[x];
		  C[j] = (C[j]>>8) | CMap[x];
#ifdef DEBUG_SCAN
		  printf(" %4d: %s %0*llx %0*llx\n",i,fmer[x],Mmer/2,N[j],Mmer/2,C[j]);
#endif
		}
	  for (i = Mmer-1; i < Kmer; i++)
		{ x = ((x<<2) | s[i]) & 0xff;
		  j = (i & 0x3);
		  n = N[j] = ((N[j]<<8) & NMask) | TMap[x];
		  c = C[j] = (C[j]>>8) | CMap[x];
#ifdef DEBUG_SCAN
		  printf(" %4d: %s %0*llx %0*llx\n",i,fmer[x],Mmer/2,n,Mmer/2,c);
#endif
		  if (n < c)
			mz = n;
		  else
			mz = c;
		  mzr[i] = mz;
		  if (mz <= mb)
			{ mi = i;
			  mb = mz;
			}
		}
	  for (i = Kmer; i < len; i++)
		{ x = ((x<<2) | s[i]) & 0xff;
		  j = (i & 0x3);
		  n = N[j] = ((N[j]<<8) & NMask) | TMap[x];
		  c = C[j] = (C[j]>>8) | CMap[x];
#ifdef DEBUG_SCAN
		  printf(" %4d: %s %0*llx %0*llx :: %4d %0*llx\n",i,fmer[x],Mmer/2,n,Mmer/2,c,mi,Mmer/2,mb);
#endif
		  if (n < c)
			mz = n;
		  else
			mz = c;
		  mzr[i&ModMask] = mz;
#ifdef DEBUG_SCAN
		  printf(" %4d: %s %0*llx %0*llx :: %4d %0*llx\n",i,fmer[x],Mmer/2,n,Mmer/2,c,mi,Mmer/2,mb);
#endif
		  if (i > mi+MaxSuper)
			{ ohang = (i-last)-1;

			  p = pad;
			  y = (mb >> p);
			  v = Count[y];
			  while (v < 0)
				{ p -= 2;
				  y  = ((mb >> p) & 0x3) - v;
				  v  = Count[y];
				}

			  transmit(v,packs,r+last,ohang); // Super-mer to i-1;
			  last = i;

			  mb = NMask;
			  j  = mi+1;
			  for (mi = -1; j <= i; j++)
				{ mz = mzr[j&ModMask];
				  if (mz <= mb)
					{ mi = j;
					  mb = mz;
					}
				}
#ifdef DEBUG_SCAN
			  printf("%*s:: %4d %0*llx  Forced\n",4*Mmer+11,"",mi,Mmer/2,mb);
#endif
			}
		  else if (mz < mb)
			{ ohang = (i-last)-1;

			  p = pad;
			  y = (mb >> p);
			  v = Count[y];
			  while (v < 0)
				{ p -= 2;
				  y  = ((mb >> p) & 0x3) - v;
				  v  = Count[y];
				}

			  transmit(v,packs,r+last,ohang); // Super-mer to i-1;
			  last = i;
			  mi = i;
			  mb = mz;
#ifdef DEBUG_SCAN
			  printf("%*s:: %4d %0*llx  New Min\n",4*Mmer+11,"",mi,Mmer/2,mb);
#endif
			}
		}
	  ohang = (len-last)-1;

	  p = pad;
	  y = (mb >> p);
	  v = Count[y];
	  while (v < 0)
		{ p -= 2;
		  y  = ((mb >> p) & 0x3) - v;
		  v  = Count[y];
		}

	  transmit(v,packs,r+last,ohang); // Super-mer to i-1;
  
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

	printf("   Transmit %d %d to bucket 1\n       ",len,ohang);
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
#ifdef DEBUG_RECIEVE
  printf("\n   Next Packet of %lld bits",totbits);
  if (bundle->erem == 64)
	printf(" covering %lld uint64 words\n\n",totbits>>6);
  else
	printf(" covering %lld uint64 words\n\n",(totbits>>6)+1);
  fflush(stdout);
#endif
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
  printf("   Transmit %d %d to bucket 1\n       ",ohang+Kmer,ohang);
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

  //  Get the hash (1st 32bp or k-mer + zero-padding if k < 32) of the cannonical k-mer at the
  //    given offset from finger.  Return +1 if the hash is from the forward direction,
  //    -1 if from the complement direction, and 0 if the two are equal.

int Get_Hash(uint64 *hash, uint64 *finger, uint64 offset)
{ uint64 *ptr;
  int     rem;
  uint64  n, c;

  rem = 64-(offset & 0x3f);
  ptr = finger + (offset>>6);
  n = get_norm64(rem,ptr);
  if (rem <= KBits)
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

  rem = 64 - (offset & 0x3f);
  ptr = finger + (offset>>6);
  if (dir <= 0)
	{ cem  = KBits-rem;
	  ctr  = ptr + (cem>>6);
	  cem  = 64 - (cem & 0x3f);
	}
  ptr += 1;

  i = 2;
  if (dir == 0)
	{ for (; i < KWords; i++)
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
	  *super = (get_comp64(cem,ctr) & KMask);
	  return (-1);
	}

  for ( ; i < KWords; i++)
	*super++ = get_norm64(rem,ptr++);
  *super = (get_norm64(rem,ptr) & KMask);
  return (1);
}
