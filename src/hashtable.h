#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <inttypes.h>
#include <Rcpp.h>

// struct ptr_size {
//   int *ptr;
//   int count;
// };
//
// typedef struct ptr_size ptr_size_t;

struct KeyValue {
  uint64_t key;
  uint64_t pos;
  int count;
};

typedef struct KeyValue keyvalue;

struct HashTable {
  int size;
  uint64_t last_pos;
  int *ids;
  keyvalue *table;
};

typedef struct HashTable hashtable;


uint64_t hashCode(uint64_t x, const int size);
hashtable *ht_create( const int size , const int idSize);
void ht_insert(hashtable *ht, const uint64_t key, const int* ids, const int count);
int *ht_get(hashtable* ht, const uint64_t key, int *count);
int *getAround(hashtable *ht, const double x, const double y, const double dist, const int mult, const int64_t maxX, const int64_t maxY, int *size);
