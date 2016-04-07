#
## a cylinder
#
algebraic3d

# cut cylinder by planes:

solid fincyl = cylinder ( 50, 0, 0; -50, 0, 0; 4)
	and plane (-48, 0, 0; -1, 0, 0)
	and plane (48, 0, 0; 1, 0, 0);

tlo fincyl;
