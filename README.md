We reduce the Ketje Sr v1 to a 200-bit state, and give an 5-round experiment using 16-dimension cube. The experiment follows the attack on 6-round initialization of Ketje Sr v1, using the same linear structure and dynamic cubes shown in Figure 12 of the paper. As the state size is reduce to 200-bit, the full key size is reduced to 64-bit key. To accelerate the test, we set the key in A[3,0] and half of A[3,1] to be zero. So totally (64-8-4=52) bits key are secret. Note that half of A[0,0] and half of A[3,1] are paddings. Without affacting the test, we simply set all paddings to be zero. 

In a PC using a single core (Intel(R)Core(TM)i7-4790 CPU@3.60GHz) with Visual studio 2010 Release x64 platform, it cost about 2 hours to recover a 24-bit equivalent key. This experiment confirms correctness our attacks.

The experiment result (The key is randomly generated.):
24-bit Equivalent Right Key (Key[0]^Key[5],Key[2]^Key[7],Key[4]): 0x21,0xd9,0x73

24-bit Equivalent Candidate Key (Key[0]^Key[5],Key[2]^Key[7],Key[4]): 0x21,0xd9,0x73


