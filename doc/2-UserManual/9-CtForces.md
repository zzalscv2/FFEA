
Constant forces {#ctforces}
===========================

Constant forces can be added onto nodes or full blobs. In order to do that,
 one needs to activate the calculation of the constant forces:

    < calc_ctforces = 1 >

in the input .ffea file. In addition to that, ` <system> ` has to include 
 the ` <interactions> ` block with the ` <ctforces> ` block, where the 
 ` ctforces_fname ` needs to be specified. Thus, a valid input ffea file 
 simulating two blobs with constant forces would have the structure:

     <param>
        ...
        <calc_springs=1>
        ...
     </param>
     <system>
        <blob>
           ...
        </blob>
        <blob>
           ...
        </blob>
        <interactions>
           <ctforces>
              <springs_fname=myspringsfile.ctforces>
           </ctforces>
        <interactions>
     </system>

Full description of the format of the input ` ctforces_fname ` file can be found 
 [here](\ref ifctforces)


