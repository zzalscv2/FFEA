#include "SparseMatrixTypes.h" 

		sparse_entry_sources::sparse_entry_sources() {
			num_sources = 0;
			sources = NULL;
		}

		sparse_entry_sources::~sparse_entry_sources() {
			num_sources = 0;
			delete[] sources;
			sources = NULL;
		}

		int sparse_entry_sources::init(int num_sources) {
			this->num_sources = num_sources;
			sources = new scalar*[num_sources];

			if(sources == NULL) {
				FFEA_ERROR_MESSG("Could not allocate memory (for sources array)\n")
			}

			return FFEA_OK;
		}

		void sparse_entry_sources::set_source(int i, scalar *s)
		{
			sources[i] = s;
		}

		scalar sparse_entry_sources::sum_all_sources()
		{
			scalar sum = 0.0;
			for(int i = 0; i < num_sources; i++) {
				sum += *(sources[i]);
			}
			return sum;
		}

		int sparse_entry_sources::get_num_sources() {
			return num_sources;
		}

