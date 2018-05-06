
#ifndef TABLE_H_
#define TABLE_H_ 1

/**
	Table object:  a table can be in three states,

	o) if t_ is a valid pointer, num_enties >=0, alloc_size > 0.
	o) t_ is NULL and both num_entries and alloc_size is 0.

	Testing for alloc_size == 0 is the same as an unitialized Table.

	Table* Table_new(int n);
	int Table_insert(Table *T, void *S);
	void*  Table_get(const Table *T, int i);
	void* Table_remove(Table *T, int i);
	int Table_num_entries(const Table *T);
	int Table_size(const Table *T);

	void Table_dump(const Table *T);

*/
typedef struct{
	void **t_;
	int num_entries_;
	int alloc_size_;
} Table;



/**
	Create a new table with a specific number of elements.

	@param n the number of items to start the list with.
	@return a new Table, or void pointer, if allocation was not successful.	
*/
Table* Table_new(int n);



/**
	Insert object into table and return table index.
	
	@return -1 if position not found.
*/	
int Table_insert(Table *T, void *S);




/**
	@param T table to look into.
	@param i location within table to retrieve item.
	@return new location if successful, -1 otherwise.
*/
void*  Table_get(const Table *T, int i);



void* Table_remove(Table *T, int i);


int Table_num_entries(const Table *T);

/**
	@return allocated size (number of entries), NOT how many are active.
*/
int Table_size(const Table *T);

/**
	Destroy table.  (Actual items are NOT destoryed, so once table is 
	destroyed, list of items are no longer available.
*/
void Table_delete(Table *T);

void Table_dump(const Table *T);

#endif


/* TABLE_H_ */
