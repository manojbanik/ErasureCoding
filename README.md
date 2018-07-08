# ErasureCoding
This program fragments a file(any type of file) into several (user can choose) pieces and then reconstruct the file even some some fragments are lost or damaged based on Reed-Solomon coding.
A block of data can be divided into n sub-blocks.  Lets call them “Data fragments” D1, D2, D3, … , Dn.     
There are some “Checksum fragments” whose contents will be calculated from the contents of the Data fragments. Checksum fragments are C1, C2, C3, … , Cm 
So we have a total of (n+m) Data & Checksum fragments.
The goal is to define the calculation of each such Checksum fragments in a way so that the loss of any m Fragments should still allow full recovery of the original Data block.
That means any 'n' out of (n+m) fragments are capable of reconstructing the original data block/file.
So the system can survive up to 'm' attacks(considering every attack damages one block).

# Have any question? Please feel free to contact/comment.
