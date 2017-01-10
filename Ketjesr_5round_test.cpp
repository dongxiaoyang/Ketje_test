#include <stdio.h>
#include <conio.h>
#include <string.h>
#include <time.h>
#include <stdlib.h>

typedef unsigned char UINT8;
typedef signed long long int INT64;

#define random(x) (rand())%x;
#define nrRounds 5
UINT8 KeccakRoundConstants[nrRounds];//these are constant,
#define nrLanes 25
unsigned int KeccakRhoOffsets[nrLanes];//these are constant,
#define nrMessage 0x8000

void KeccakPermutationOnWords(UINT8 *state);
void theta(UINT8 *A);
void rho(UINT8 *A);
void pi(UINT8 *A);
void chi(UINT8 *A);
void iota(UINT8 *A, unsigned int indexRound);



void KeccakPermutationOnWords(UINT8 *state)
{
    unsigned int i;

    for(i=0; i<nrRounds; i++) {
        theta(state);
        rho(state);
        pi(state);
        chi(state);
        iota(state, i);
    }
}

void KeccakPermutationOnWords_withkey(UINT8 *state,UINT8 *key)
{
    unsigned int i;

	state[0]=key[0];
	state[1]=key[1];
	state[2]=key[2];
	state[3]=key[3];
	state[4]=key[4];
	state[5]=key[5];
	state[6]=key[6];
	state[7]=key[7];
	state[8]=key[8];
    for(i=0; i<nrRounds; i++) {
        theta(state);
        rho(state);
        pi(state);
        chi(state);
        iota(state, i);
    }
}


// used in theta
#define index(x, y) (((x)%5)+5*((y)%5))
//used in theta
#define ROL32(a, offset) ((offset != 0) ? ((((UINT8)a) << offset) ^ (((UINT8)a) >> (8-offset))) : a)

void theta(UINT8 *A)
{
    unsigned int x, y;
    UINT8 C[5], D[5];//C are the Xors of the five bits in every column. D are the Xors of the ten bits in right-behind column and right column

    for(x=0; x<5; x++) {
        C[x] = 0;
        for(y=0; y<5; y++)
            C[x] ^= A[index(x, y)];
    }
    for(x=0; x<5; x++)
        D[x] = ROL32(C[(x+1)%5], 1) ^ C[(x+4)%5];
    for(x=0; x<5; x++)
        for(y=0; y<5; y++)
            A[index(x, y)] ^= D[x];
}

void rho(UINT8 *A)
{
    unsigned int x, y;

    for(x=0; x<5; x++) for(y=0; y<5; y++)
        A[index(x, y)] = ROL32(A[index(x, y)], KeccakRhoOffsets[index(x, y)]);
}

void pi(UINT8 *A)
{
    unsigned int x, y;
    UINT8 tempA[25];

    for(x=0; x<5; x++) for(y=0; y<5; y++)
        tempA[index(x, y)] = A[index(x, y)];
    for(x=0; x<5; x++) for(y=0; y<5; y++)
        A[index(0*x+1*y, 2*x+3*y)] = tempA[index(x, y)];//learn from this!
}

void chi(UINT8 *A)
{
    unsigned int x, y;
    UINT8 C[5];

    for(y=0; y<5; y++) {
        for(x=0; x<5; x++)
            C[x] = A[index(x, y)] ^ ((~A[index(x+1, y)]) & A[index(x+2, y)]);
        for(x=0; x<5; x++)
            A[index(x, y)] = C[x];
    }
}

void iota(UINT8 *A, unsigned int indexRound)
{
    A[index(0, 0)] ^= KeccakRoundConstants[indexRound];
}

int LFSR86540(UINT8 *LFSR)
{
    int result = ((*LFSR) & 0x01) != 0;
    if (((*LFSR) & 0x80) != 0)
        // Primitive polynomial over GF(2): x^8+x^6+x^5+x^4+1
        (*LFSR) = ((*LFSR) << 1) ^ 0x71;
    else
        (*LFSR) <<= 1;
    return result;
}

void KeccakInitializeRoundConstants()
{
    UINT8 LFSRstate = 0x01;
    unsigned int i, j, bitPosition;

    for(i=0; i<nrRounds; i++) {
        KeccakRoundConstants[i] = 0;
        for(j=0; j<7; j++) {
            bitPosition = (1<<j)-1; //2^j-1
            if (LFSR86540(&LFSRstate))
                KeccakRoundConstants[i] ^= (UINT8)1<<bitPosition;
        }
    }
}

void KeccakInitializeRhoOffsets()
{
    unsigned int x, y, t, newX, newY;

    KeccakRhoOffsets[index(0, 0)] = 0;
    x = 1;
    y = 0;
    for(t=0; t<24; t++) {
        KeccakRhoOffsets[index(x, y)] = ((t+1)*(t+2)/2) % 8;
        newX = (0*x+1*y) % 5;
        newY = (2*x+3*y) % 5;
        x = newX;
        y = newY;
    }
}

void KeccakInitialize()
{
    KeccakInitializeRoundConstants();
    KeccakInitializeRhoOffsets();
}



int main(int argc,char *argv[])
{
	clock_t start,finish;
	start = clock();
	
	FILE *f;
	f=fopen("table.txt","w+b");

    UINT8 InitialState[25]={0};
	UINT8 TempState[25]={0};
	UINT8 FinalState[4]={0};
	UINT8 KeyCubesum[256][8]={0};
	UINT8 GuessCubesum[8]={0};
	UINT8 Key[9]={0};
	INT64 i,j,k,temp,temp1,counter;

	KeccakInitialize();

	unsigned int indexk[2][20];
	indexk[0][0]=12;indexk[1][0]=0;  
	indexk[0][1]=12;indexk[1][1]=1;  
	indexk[0][2]=12;indexk[1][2]=2;  
	indexk[0][3]=12;indexk[1][3]=3; 
	indexk[0][4]=12;indexk[1][4]=4; 
	indexk[0][5]=12;indexk[1][5]=5; 
	indexk[0][6]=12;indexk[1][6]=6; 
	indexk[0][7]=12;indexk[1][7]=7; 


	unsigned int index[2][16];
	index[0][0]=14;index[1][0]=0;  
	index[0][1]=14;index[1][1]=1;  
	index[0][2]=14;index[1][2]=2;  
	index[0][3]=14;index[1][3]=3; 
	index[0][4]=14;index[1][4]=4; 
	index[0][5]=14;index[1][5]=5; 
	index[0][6]=14;index[1][6]=6; 
	index[0][7]=14;index[1][7]=7; 
	index[0][8]=19;index[1][8]=0; 
	index[0][9]=19;index[1][9]=1;  
	index[0][10]=19;index[1][10]=2;  
	index[0][11]=19;index[1][11]=3; 
	index[0][12]=19;index[1][12]=4; 
	index[0][13]=19;index[1][13]=5; 
	index[0][14]=19;index[1][14]=6; 
	index[0][15]=19;index[1][15]=7;
	//preprocess
	for(i=0;i<(1<<8);i++)
	{
		for(j=0;j<25;j++){
			InitialState[j]=0;
		}
		for(j=0;j<4;j++){
			FinalState[j]=0;
		}			
		fprintf(f,"key bits value:");
		InitialState[3] = 0;
		InitialState[4] = i;
		InitialState[8] = 0;
		fprintf(f,"%02x.",i);
		

		for(j=0;j<(1<<16);j++)//initialize tempstate
		{
			for(k=0;k<25;k++)//fresh the tempstate for the key and initial value
			{
				TempState[k]=InitialState[k];
			}
			//give the value for cube variables
			TempState[15] = (j&0xff);
			TempState[22] = ((j>>8)&0xff);
			//give the value for dynamic value
			TempState[10]=TempState[15];
			TempState[17]=TempState[22];
			KeccakPermutationOnWords((UINT8 *)TempState);
			for(k=0;k<4;k++){
				FinalState[k]^=TempState[k];
			}
		}
		for(j=0;j<4;j++){
			fprintf(f,"%02x ",FinalState[j]);
			KeyCubesum[i][j]=FinalState[j];
		}
		for(j=0;j<4;j++){
			FinalState[j]=0;
		}	
		InitialState[12]^=0xff;
		for(j=0;j<(1<<16);j++)//initialize tempstate
		{
			for(k=0;k<25;k++)//fresh the tempstate for the key and initial value
			{
				TempState[k]=InitialState[k];
			}
			//give the value for cube variables
			TempState[15] = (j&0xff);
			TempState[22] = ((j>>8)&0xff);
			//give the value for dynamic value
			TempState[10]=TempState[15];
			TempState[17]=TempState[22];
			KeccakPermutationOnWords((UINT8 *)TempState);
			for(k=0;k<4;k++){
				FinalState[k]^=TempState[k];
			}
		}
		for(j=0;j<4;j++){
			fprintf(f,"%02x ",FinalState[j]);
			KeyCubesum[i][4+j]=FinalState[j];
		}
		fprintf(f,"\n");
	}
	start = clock();
	//online
	srand(time(NULL));
	//generate a random 64-bit key for Key[9],
	for(i=0;i<9;i++)
	{
		for(j=0;j<8;j++)
		{
			temp=random(2);
			if(temp)
			{
				Key[i]^=(1<<j);
			}
		}
	}
	//note that we add padding=0000 in MSB of Key[0] and LSB of Key[8].
	Key[0]&=0x0f;
	Key[3]=0;
	Key[8]=0;
	
	printf("24-bit Equivalent Right Key (Key[0]^Key[5],Key[2]^Key[7],Key[4]): 0x%02x,0x%02x,0x%02x\n",Key[0]^Key[5],Key[2]^Key[7],Key[4]);
	//key guessing phase 
	
	for(i=0;i<(1<<16);i++){
		//if(((i)&0xff)!=(Key[0]^Key[5]))continue;
		//if(((i>>8)&0xff)!=(Key[2]^Key[7]))continue;
		for(j=0;j<25;j++){
			InitialState[j]=0;
		}
		for(j=0;j<4;j++){
			FinalState[j]=0;
		}			
	
		for(j=0;j<(1<<16);j++)//initialize tempstate
		{
			for(k=0;k<25;k++)//fresh the tempstate for the key and initial value
			{
				TempState[k]=InitialState[k];
			}
			//give the value for cube variables
			TempState[15] = (j&0xff);
			TempState[22] = ((j>>8)&0xff);
			//give the value for dynamic value
			TempState[10]=TempState[15]^(i&0xff);
			TempState[17]=TempState[22]^((i>>8)&0xff);
			KeccakPermutationOnWords_withkey((UINT8 *)TempState,Key);
			for(k=0;k<4;k++){
				FinalState[k]^=TempState[k];
			}
		}
		for(j=0;j<4;j++){
			GuessCubesum[j]=FinalState[j];
		}
		for(j=0;j<4;j++){
			FinalState[j]=0;
		}	
		InitialState[12]^=0xff;
		for(j=0;j<(1<<16);j++)//initialize tempstate
		{
			for(k=0;k<25;k++)//fresh the tempstate for the key and initial value
			{
				TempState[k]=InitialState[k];
			}
			//give the value for cube variables
			TempState[15] = (j&0xff);
			TempState[22] = ((j>>8)&0xff);
			//give the value for dynamic value
			TempState[10]=TempState[15]^(i&0xff);
			TempState[17]=TempState[22]^((i>>8)&0xff);
			KeccakPermutationOnWords_withkey((UINT8 *)TempState,Key);
			for(k=0;k<4;k++){
				FinalState[k]^=TempState[k];
			}
		}
		for(j=0;j<4;j++){
			GuessCubesum[4+j]=FinalState[j];
		}
		for(j=0;j<256;j++){
			if(memcmp(&KeyCubesum[j][0],&GuessCubesum[0],8*sizeof(UINT8))==0){
				for(int bb=0;bb<8;bb++)
				printf("%02x ",KeyCubesum[j][bb]);
				printf(" = ");
				for(int bb=0;bb<8;bb++)
				printf("%02x ",GuessCubesum[bb]);
				printf("\n");
				
				printf("24-bit Equivalent Candidate Key (Key[0]^Key[5],Key[2]^Key[7],Key[4]): 0x%02x,0x%02x,0x%02x\n",(i&0xff),((i>>8)&0xff),j);
	
			}
		}
	}


	fclose(f);
	finish = clock();
	printf("time=%f\n",(double)(finish-start)/CLK_TCK);
	printf("Done!\n");
    getch();
	
    return 0;
}


