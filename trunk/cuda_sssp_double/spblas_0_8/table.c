#include <stdlib.h>
#include <stdio.h>
#include "table.h"

#ifndef NULL
#define NULL 0
#endif

/*
#ifdef DEBUG
#define DEBUG_OP(x) x
#else
#define DEBUG_OP(x)
#endif
*/

#define DEBUG_PRINT 0

Table* Table_new(int num)
{
	Table *T = NULL;

	if (num < 1)
		return NULL;

	T = (Table *) malloc(sizeof(Table));
	if (T != NULL)
	{
		T->num_entries_ = 0;
		T->alloc_size_ =0;
		T->t_ = (void **) malloc(sizeof(void *) * num);
		if (T->t_ != NULL)
		{
			int i=0;
			T->alloc_size_ = num;
			for (; i<T->alloc_size_; i++)
				T->t_[i] = NULL;
		}
	}
	return T;
}

void* Table_remove(Table *T, int i)
{
	void *temp = NULL;

	if (T == NULL) return NULL;
	if (i < 0 ) return NULL;
	if (i >= T->alloc_size_) return NULL;

	temp = T->t_[i];
	T->t_[i] = NULL;
	T->num_entries_--;
	return temp;
}

int Table_size(const Table *T)
{
	return (T != 0) ? T->alloc_size_ : 0;

}

int Table_num_entries(const Table *T)
{
	return (T != 0) ? T->num_entries_ : 0;

}

void Table_delete(Table *T)
{
	if (T==NULL)
		return;

	free(T->t_);
	free(T);
}



void* Table_get(const Table *T, int i)
{
	if (T==NULL) return NULL;
	if (i < 0 ) return NULL;
	if (i >= T->alloc_size_) return NULL;

	return T->t_[i];
}

static int Table_grow_(Table *T, int newsize)
{
	int i = 0;
	void **bt = (void **) realloc(T->t_, newsize*sizeof(void*));

	if (bt == NULL)
		return -1;

	T->t_ = bt;
	/* initialize new enties */
	for (i=T->alloc_size_; i< newsize; i++)
		T->t_[i] = NULL;
	T->alloc_size_ = newsize;

	return 0;
}

int Table_insert(Table *T, void *S)
{
	int i=0;			/* counter used for searching empty slot. */

	/*
     * First, find out if we need to grow current table.
	 * If the table is empty, initialize it with one available slot.
	 * Otherwise, if the table is already full, we will need to
	 * grow it.  (We double the size of the table.)
	 *
	 * Now we can find an available slot by looking for
	 * NULL entires in the table.  (The grow() and init()
	 * routines both set unused entries in table to NULL.)
	 *
	 * If for some reason, we still can't find available
	 * slots in the table, then the table has somehow been
	 * corrupted.  Best thing we can do at this point, is
	 * to flag this situtation, and return without modification
	 * to input paramters.
	 *
	*/


	/* premature return if table is corrupt */
	if (T == NULL || 
		T->t_ == NULL || 
		T->alloc_size_ == 0)
			return -1;

	if (DEBUG_PRINT)
		printf("Table_insert: passed NULL table test.\n");

	/* if table is full, expand it. */
	if (!(T->num_entries_ < T->alloc_size_))
	{

		if (Table_grow_(T, 2*T->alloc_size_) != 0)
			return -1;
	}

	if (DEBUG_PRINT)
		printf("Table_insert: passed grow table test.\n");


   /*
    * At this point, there is an available slot -- we just need to find it.
    *
    * Performance hack: don't start looking at the begining of the list.
    * This optimizes for the common case where one is adding and removing
    * matrices at the end of list (i.e. using the table as a stack).
    * Otherwise it becomes an O(n), rather than an O(1) insertion cost.
   */


	for (i=0; i< T->alloc_size_; i++)
	{
		int pos = (T->num_entries_ + i ) % T->alloc_size_;
		if (T->t_[pos] == NULL)
		{
			/* found an empty (available) slot. */
			if (DEBUG_PRINT)
				printf("Table_insert: found empty slot [%d].\n", pos);
			T->t_[pos] = S;
			T->num_entries_++;
			return pos;
		}
	
	}


	/*  if no position found, we have a corrupt table (or ran out of memory.)
	 *  best thing to do is to return an error condition.
	*/
	return -1;

}

void Table_dump(const  Table *T)
{
	int i;

	printf("T: %d\n", (int) T);
	if (T==NULL) return;

	printf("Alloc size:  %d\n", T->alloc_size_);
	printf("Num entries: %d\n", T->num_entries_);
	for (i=0; i<T->alloc_size_; i++)
	{
		if (T->t_[i] == NULL)
			printf("[%d]  NULL \n",i);
		else
			printf("[%d]   %d \n", i, (int) T->t_[i]);
	}
	printf("\n");
}

	


