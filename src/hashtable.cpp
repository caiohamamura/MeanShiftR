#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <inttypes.h>
#include <Rcpp.h>
#include "hashtable.h"




hashtable *ht_create( const int size , const int idSize) {
  hashtable* ht = NULL;

  if( size < 1 ) return ht;
  if( ( ht = (hashtable*)malloc( sizeof( hashtable ) ) ) == NULL ) {
    return NULL;
  }

  ht->table = (keyvalue*)malloc( sizeof( keyvalue ) * size );
  ht->ids = (int*)malloc( sizeof( int ) * idSize );
  ht->size = size;
  ht->last_pos = 0;

  for (int i = 0; i < size; i++) {
    ht->table[i].count = -1;
  }

  return ht;
}

uint64_t hashCode(uint64_t x, const int size) {
  x = (x ^ (x >> 30)) * UINT64_C(0xbf58476d1ce4e5b9);
  x = (x ^ (x >> 27)) * UINT64_C(0x94d049bb133111eb);
  x = x ^ (x >> 31);
  return (x % size);
}

void ht_insert(hashtable *ht, const uint64_t key, const int* ids, const int count) {
  //get the hash
  uint64_t hashIndex = hashCode(key, ht->size);

  //move in array until an empty or deleted cell
  while(ht->table[hashIndex].count != -1) {
    if (ht->table[hashIndex].key == key) {
      break;
    }
    //go to next cell
    ++hashIndex;

    //wrap around the table
    hashIndex %= ht->size;
  }

  ht->table[hashIndex].key = key;
  ht->table[hashIndex].count = count;
  ht->table[hashIndex].pos = ht->last_pos;
  // Rcpp::Rcout <<
  //   "I inserted key: " << ht->table[hashIndex].key << std::endl;

  for (int i = 0; i < count; i++) {
    ht->ids[ht->last_pos++] = ids[i];
  }
}

int *ht_get(hashtable* ht, const uint64_t key, int *count) {
  //get the hash
  uint64_t hashIndex = hashCode(key, ht->size);

  //move in array until an empty
  while(ht->table[hashIndex].count != -1) {

    if(ht->table[hashIndex].key == key) {
      *count = ht->table[hashIndex].count;
      return &(ht->ids[ht->table[hashIndex].pos]);
    }

    //go to next cell
    ++hashIndex;

    //wrap around the table
    hashIndex %= ht->size;
  }

  return NULL;
}


int *getAround(hashtable *ht, const double x, const double y, const double dist, const int mult, const int64_t maxX, const int64_t maxY, int *oSize) {
  int n = 0;
  int k = 0;
  int64_t xMin = (x - dist);
  int64_t xMax = (x + dist);
  int64_t yMin = (y - dist);
  int64_t yMax = (y + dist);
  if (xMin < 0) xMin = 0;
  if (yMin < 0) yMin = 0;
  if (xMax > maxX) xMax = maxX;
  if (yMax > maxY) yMax = maxY;
  int *out = (int*)malloc(sizeof(int) * 500);
  for (int64_t curY = yMin; curY <= yMax; curY++) {
    for (int64_t curX = xMin; curX <= xMax; curX++) {
      uint64_t idx = (curX * mult) + curY;
      // Rcpp::Rcout << "\ridx: " << idx << "    ";
      int size = -1;
      // Rcpp::Rcout << "Getting from ht...\n";
      int* ids = ht_get(ht, idx, &size);
      n += size;
      // Rcpp::Rcout << "Reallocating to size:" << size << "...\n";
      // Rcpp::Rcout << "Putting values...\n";
      for (int it = 0; it < size;) {
        out[k++] = ids[it++];
      }
    }
  }
  *oSize = n;
  return out;
}
