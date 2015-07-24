#ifndef NODEMAP_H_INCLUDED
#define NODEMAP_H_INCLUDED

#include ""
class NodeMap {

	public:
		
		NodeMap(); 

		~NodeMap(); 
	
		int load_map_from_file(int from_conformation, to_conformation, fname);


	// Variables

	// conformation mapped from index
	int from_conformation;

	// conformation mapping to index
	int to_conformation;
};
#endif
