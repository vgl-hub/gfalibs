#ifndef MIN_SCAN

#define MIN_SCAN

#include <pthread.h>

typedef unsigned char      uint8;
typedef unsigned long long uint64;
typedef signed   long long int64;

typedef void               Distribution_Bundle;
typedef void               Scan_Bundle;

// DISTRIBUTION PHASE:

  //  To start the distribution phase, tell the package the kmer size and the number of files
  //    the supermers should be distributed over.  Also give it an initial cache of data of
  //    length clen to train on.  Based on this segment the routine will select a minimizer
  //    size to use and develop a mapping scheme, as per the one used in FastK.  The routine
  //    returns 0 if it could not produce a scheme (should never happen basically) and otherwise
  //    the size of the minimizer it will use.

int Train_Genes_Package(int kmer, int fnum, char *cache, int clen);

  //  For each thread working on a segment of the data, first call Begin_Distribution
  //    with an array of 'fnum' (see training above) system file descriptors, one for
  //    each bucket.  Then call Distribute_Sequence on successive blocks of the segment --
  //    every sequence is assumed to be terminated with an 'n' implying the last character
  //    of reads is expected to be 'n'.  When the segment has been processed call
  //    End_Distribution to flush the write buffers and free the worrking storage allocated
  //    at the start.

Distribution_Bundle *Begin_Distribution(int *fids);

void Distribute_Sequence(char *reads, int rlen, Distribution_Bundle *bundle);

void End_Distribution(Distribution_Bundle *bundle);


// COUNTING PHASE:

  //  Before entering threaded code, initialize the package with this call
  //    giving the kmer and mmer size for the run.  Both these parameters should
  //    be those used in the distribution phase.

void Init_Genes_Package(int kmer, int mmer);

  //  Begin scanning a data segment (that must be some number of complete packets)
  //    starting at data and covering size *uint64's*.  Each call to Next_Supermer delivers
  //    the next supermer in the data segement in the use-supplied array super, returning
  //    the length of the supermer.  0 is returned if there are no more supermers to process.
  //    A uint64 array of the appropriate size to hold any supermer is returned by
  //    New_Supermer_Buffer.  The 1st base pair of a supermer is in the highest 2 bits of
  //    the first uint64, and the last uint64 is 0-padded if not completely filled.  When
  //    the scan is complete, End_Supermer_Scan should be called with the scan bundle in
  //    order to free it.

Scan_Bundle *Begin_Supermer_Scan(uint64 *data, uint64 size);

uint64 *New_Supermer_Buffer();

int Next_Supermer(uint64 *supermer, Scan_Bundle *bundle);

void End_Supermer_Scan(Scan_Bundle *bundle);

  //  Given an intialized scan bundle at any point in a data segment, advance to the next packet
  //    and return a pointer to the start of said packet.  If there is no next packet then return
  //    NULL.

uint64 *Skip_To_Next_Packet(Scan_Bundle *bundle);

  //  Assumes the bundle points at the start of the overhang length of a supermer record.
  //  Get and scan past the overhang length and return the # of k-mers in the supermer,
  //  i.e. over hang length + 1;

int    Get_Kmer_Count(Scan_Bundle *bundle);

  //  Assumes the bundle points at the start of a 2-bit supermer containing count k-mers.
  //  Advance the bundle position over the supermer to the start of the next supermer record.

void   Skip_Kmers(int count, Scan_Bundle *bundle);

  //  Return an integer offset representing the precise bit position of the current bundle relative
  //    to a pointer finger into the data array.  You can get an offset before processing a supermer
  //    (but after getting its count with Get_Kmer_Count) and know that the successive k-mers are
  //    at offset+0, +2, +4, +6 ...

uint64 Current_Offset(uint64 *finger, Scan_Bundle *bundle);

  //  Get the hash (1st 32bp of k-mer + zero-padding if k < 32) of the cannonical k-mer at the
  //    given offset from finger.  Return +1 if the hash is from the forward direction,
  //    -1 if from the complement direction, and 0 if the two are equal.

int Get_Hash(uint64 *hash, uint64 *finger, uint64 offset);

  //  Get the canonical k-mer at the offset from finger whose hash is hash that comes from the end
  //    specified by dir (both obtained from Get_Hash).  1 is returned if the k-mer is in the
  //    forward direction and -1 if in the complement direction.

int Get_Canonical_Kmer(uint64 *super, int dir, uint64 hash, uint64 *finger, uint64 offset);

#endif
