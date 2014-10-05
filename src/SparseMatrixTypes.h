#ifndef SPARSEMATRIXTYPES_H_INCLUDED
#define SPARSEMATRIXTYPES_H_INCLUDED

typedef struct
{
	/* The column index of this entry in the original matrix */
	int column_index;

	/* The value of this entry */
	scalar val;
} sparse_entry;

class sparse_entry_sources
{
	public:
		sparse_entry_sources() {
			num_sources = 0;
			sources = NULL;
		}

		~sparse_entry_sources() {
			num_sources = 0;
			delete[] sources;
			sources = NULL;
		}

		int init(int num_sources) {
			this->num_sources = num_sources;
			sources = new scalar*[num_sources];

			if(sources == NULL) {
				FFEA_ERROR_MESSG("Could not allocate memory (for sources array)\n")
			}

			return FFEA_OK;
		}

		void set_source(int i, scalar *s)
		{
			sources[i] = s;
		}

		scalar sum_all_sources()
		{
			scalar sum = 0.0;
			for(int i = 0; i < num_sources; i++) {
				sum += *(sources[i]);
			}
			return sum;
		}

		int get_num_sources() {
			return num_sources;
		}

	private:
		int num_sources;
		scalar **sources;
};

#endif
