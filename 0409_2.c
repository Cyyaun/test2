#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h> 
#define PI 3.1415926

int row_num = 80;

typedef struct AtomRecord {
	struct AtomRecord *anext;
	int AtomNumber;
	char AtomName[4];
	char Alternate;
	char ResName[4];
	char Chain;
	int ResNumber;
	double Coordx,Coordy,Coordz;
	float Occupancy;
	float Bfactor;
}Atom;

typedef struct ResidueRecord {
	struct ResidueRecord *rnext;
	Atom *Natom;
	Atom *CAatom;
	Atom *Catom;
	Atom *Oatom;
	double TCO, Kappa, Alpha, phi, psi;
	Atom *CBatom;// 
	char ResName3[4];// 
	char ResName1;//
	int ResType;//
	int ResNumber;// 
	int AtomCount;//
}Residue;

/*
typedef struct PeptideRecord {
	struct PeptideRecord *pnext;
	Atom *FirstAtom;
	Atom *LastAtom;
	Residue *FirstResidue;
	Residue *LastResidue;
	//char Chain;
	//int AtomCount;
	//int ResCount;
} Peptide;
*/

typedef struct ProteinRecord {
	Atom *FirstAtom;
	Atom *LastAtom;
	Residue *FirstResidue;
	Residue *LastResidue;
	//Peptide *FirstPeptide;
	//Peptide *LastPeptide;
	char PDBname[5];//
	int AtomCount;//
	int ResCount;//
	int PepCount;//
} Protein;

Protein TestProtein;

float Distance(float x1,float y1,float z1,float x2,float y2,float z2)
{
	float x,y,z,d;
	
	x=x2-x1;
	y=y2-y1;
	z=z2-z1;
	d=sqrt(x*x+y*y+z*z);
	
	return(d); 
}

void DistanceMatrix(int len)
{
	FILE *fptr3;
	fptr3 = fopen("./dm.txt","w");
	if(fptr3 == NULL)
	{
		printf ("Can't write to the file path: ./dm.txt \n"); 
		exit(0);
	}
	
	int temp = len*len;
	float disarray[temp];
	
	Residue *TempResidue1,*TempResidue2;
	int i=0,j=0;
	double dis;
	
	for(TempResidue1=TestProtein.FirstResidue;TempResidue1->rnext;TempResidue1=TempResidue1->rnext)
	{
		j=0;
		for(TempResidue2=TestProtein.FirstResidue;TempResidue2->rnext;TempResidue2=TempResidue2->rnext)
		{
			dis = Distance(\
			 TempResidue1->Catom->Coordx,TempResidue1->Catom->Coordy,TempResidue1->Catom->Coordz,\
			 TempResidue2->Catom->Coordx,TempResidue2->Catom->Coordy,TempResidue2->Catom->Coordz);		
			disarray[i*len+j] = round(dis);
			j+=1;
		}
		i+=1;
	}
	
	fprintf(fptr3,"   ");
	for(i=0;i<len;i++)
	{
		if((i+1)%10==0)
			fprintf(fptr3,"%c",'A'+(i+1)/10-1);		
		else
			fprintf(fptr3,"%d",(i+1)%10);
		
	}
	fprintf(fptr3,"\n");
	
	for(i=0;i<len;i++)
	{
		fprintf(fptr3,"%2d ",i+1);
		for(j=0;j<len;j++)
		{
			if(disarray[i*len+j]>9)
				fprintf(fptr3,".");
			else
				fprintf(fptr3,"%.0f",disarray[i*len+j]);
		}
		fprintf(fptr3,"\n");
	}
	
	fclose(fptr3);
}


float cosAngle(float x1,float y1,float z1,float x2,float y2,float z2,float x3,float y3,float z3,float x4,float y4,float z4)
{
	float r1x,r1y,r1z;
	float r2x,r2y,r2z;
	float bi,bk;
	float cosangle;
	
	r1x=x1-x2;
	r1y=y1-y2;
	r1z=z1-z2;

	r2x=x3-x4;
	r2y=y3-y4;
	r2z=z3-z4;	
	
	bi=Distance(r1x,r1y,r1z,0,0,0);
	bk=Distance(r2x,r2y,r2z,0,0,0);
	
	if((bi==0)||(bk==0))
		cosangle=1;
	else
		cosangle=(r1x*r2x+r1y*r2y+r1z*r2z)/(bi*bk);
	
	return(cosangle);
}

float Angle(float x1,float y1,float z1,float x2,float y2,float z2,float x3,float y3,float z3,float x4,float y4,float z4)
{
	float r1x,r1y,r1z;
	float r2x,r2y,r2z;
	float bi,bk;
	float cosangle,ang;
	
	r1x=x1-x2;
	r1y=y1-y2;
	r1z=z1-z2;

	r2x=x3-x4;
	r2y=y3-y4;
	r2z=z3-z4;	
	
	bi=Distance(r1x,r1y,r1z,0,0,0);
	bk=Distance(r2x,r2y,r2z,0,0,0);
	
	if((bi==0)||(bk==0))
		cosangle=1;
	else
		cosangle=(r1x*r2x+r1y*r2y+r1z*r2z)/(bi*bk);
	
	//acos給cosangle回傳θ(徑度π)	
	ang = (acos(cosangle)*180)/PI;
	
	return(ang);
}

float TorsionAngle(float x1,float y1,float z1,float x2,float y2,float z2,float x3,float y3,float z3,float x4,float y4,float z4)
{
	float r1x,r1y,r1z;
	float r2x,r2y,r2z;
	float p1x,p1y,p1z;
	float p2x,p2y,p2z;
	float q1x,q1y,q1z;
	float cosang,ang,sign;
	
	r1x=x1-x2;
	r1y=y1-y2;
	r1z=z1-z2;
	r2x=x3-x2;
	r2y=y3-y2;
	r2z=z3-z2;
	p1x=r1y*r2z-r1z*r2y;
	p1y=r1z*r2x-r1x*r2z;
	p1z=r1x*r2y-r1y*r2x;
	
	r1x=x2-x3;
	r1y=y2-y3;
	r1z=z2-z3;
	r2x=x4-x3;
	r2y=y4-y3;
	r2z=z4-z3;
	p2x=r1y*r2z-r1z*r2y;
	p2y=r1z*r2x-r1x*r2z;
	p2z=r1x*r2y-r1y*r2x;
	
	q1x=p1y*p2z-p1z*p2y;
	q1y=p1z*p2x-p1x*p2z;
	q1z=p1x*p2y-p1y*p2x;
	r2x=x3-x2;
	r2y=y3-y2;
	r2z=z3-z2;
	
	cosang=cosAngle(p1x,p1y,p1z,0,0,0,p2x,p2y,p2z,0,0,0);
	ang = (acos(cosang)*180)/PI;
	sign = q1x*r2x+q1y*r2y+q1z*r2z;
	if(sign>=0)
		return(ang);
	else
		return(-ang);

}

void ReadPDB(char *PDBID)
{
	FILE *fptr;
	char filepath[100], buffer[100], section[6];
	char atomn[6], resn[5], coord[9], occ[7], bf[7];
	Atom *NewAtom;
	Atom *NAtom=NULL, *CaAtom=NULL, *CAtom=NULL, *OAtom=NULL;
	Residue *NewResidue;
	
	sprintf(filepath,"./%s.pdb",PDBID);
	fptr = fopen(filepath,"r");
	if(fptr == NULL)
	{
		printf ("Can't open the file path: %s\n",filepath); 
		exit(0);
	}
	
	int resnum=0;
	while(fgets(buffer,row_num+2,fptr) != NULL) //包含\0+\n 
	{
		printf("%s\n",buffer);
		//buffer從0開始共6字元複製到section中 
		strncpy(section,&buffer[0],6);
		section[6]='\0';
			
		//檢查section是否為"ATOM  "，是則回傳0 
		if(strcmp(section,"ATOM  ")==0)
		{ 
			NewAtom = (Atom*)malloc(sizeof(Atom));	//新list
			NewAtom->anext = NULL;
			
			//AtomNumber:7-11
			strncpy(atomn,&buffer[6],5);
			atomn[5]='\0';	
			NewAtom->AtomNumber = atoi(atomn);
			
			//AtomName:13-16
			strncpy(NewAtom->AtomName,&buffer[13],3);	
			NewAtom->AtomName[3]='\0';
			
			//Alternate:17
			NewAtom->Alternate = buffer[16];
			
			//ResName:18-20
			strncpy(NewAtom->ResName,&buffer[17],3);	
			NewAtom->ResName[3]='\0';
			
			//Chain:22
			NewAtom->Chain = buffer[21];
			
			//ResNumber:23-26
			strncpy(resn,&buffer[22],4);	
			NewAtom->ResNumber = atoi(resn);
			
			//Coordx:31-38,39-46,47-54
			strncpy(coord,&buffer[30],8);
			coord[8]='\0';	
			NewAtom->Coordx = atof(coord);
			strncpy(coord,&buffer[38],8);
			coord[8]='\0';	
			NewAtom->Coordy = atof(coord);
			strncpy(coord,&buffer[46],8);
			coord[8]='\0';	
			NewAtom->Coordz = atof(coord);
			
			//Occupancy:55-60
			strncpy(occ,&buffer[54],6);
			occ[6]='\0';	
			NewAtom->Coordz = atof(occ);
			
			//Bfactor:61-66
			strncpy(bf,&buffer[60],6);
			bf[6]='\0';	
			NewAtom->Coordz = atof(bf);
			
			
			//if atom is the first atom
			if(TestProtein.LastAtom == NULL)
			{
				TestProtein.FirstAtom = NewAtom;
				TestProtein.LastAtom = NewAtom;
			}
			//if not first then link
			else
			{
				TestProtein.LastAtom->anext = NewAtom;
				TestProtein.LastAtom = NewAtom;
			}
			
			//Atom is N 
			if((NewAtom->AtomName[0]=='N')&&(NewAtom->AtomName[1]==' '))
			{
				NAtom = NewAtom;
			}
			//Atom is CA 
			else if((NewAtom->AtomName[0]=='C')&&(NewAtom->AtomName[1]=='A'))
			{
				CaAtom = NewAtom;
			}
			//Atom is C 
			else if((NewAtom->AtomName[0]=='C')&&(NewAtom->AtomName[1]==' '))
			{
				CAtom = NewAtom;
			}
			//Atom is O 
			else if((NewAtom->AtomName[0]=='O')&&(NewAtom->AtomName[1]==' '))
			{
				OAtom = NewAtom;
			}
			
			//四種Atom's ResNumber皆一樣且不為NULL，則四種Atom屬於同一個residue且紀錄完整 
			
			if((NAtom!=NULL)&&(CaAtom!=NULL)&&(CAtom!=NULL)&&(OAtom!=NULL)\
			&&(NAtom->ResNumber==CaAtom->ResNumber)\
			&&(NAtom->ResNumber==CAtom->ResNumber)\
			&&(NAtom->ResNumber==OAtom->ResNumber)
			&&(NAtom->ResNumber!=resnum))\
			{
				NewResidue = (Residue*)malloc(sizeof(Residue));	//新list
				NewResidue->Natom = NAtom;
				NewResidue->CAatom = CaAtom;
				NewResidue->Catom = CAtom;
				NewResidue->Oatom = OAtom;
				NewResidue->rnext = NULL;
				
				printf("---ResNumber:%d---\n ",OAtom->ResNumber);
				
				//if Residue is the first Residue
				if(TestProtein.LastResidue == NULL)
				{
					TestProtein.FirstResidue = NewResidue;
					TestProtein.LastResidue = NewResidue;
				}
				//if not first then link
				else
				{
					TestProtein.LastResidue->rnext = NewResidue;
					TestProtein.LastResidue = NewResidue;
					resnum+=1; 
				} 
			}
			
		}
	}
	printf("Residue count:%d\n",resnum);
	fclose(fptr);	
}

void CalTCO()
{
	Residue *TempResidue;
	
	//前一個residue賦值
	TestProtein.FirstResidue->TCO = 0;
	
	//from first residue move to last residue
	for(TempResidue=TestProtein.FirstResidue;TempResidue->rnext;TempResidue=TempResidue->rnext)
	{
		//每一個胺基酸都要跟前一個比較(i-1 vs i)，但最後存取TempResidue->rnext
		// i = TempResidue->rnext ; i-1 = TempResidue ;
		TempResidue->rnext->TCO = cosAngle\
		(TempResidue->rnext->Catom->Coordx,TempResidue->rnext->Catom->Coordy,TempResidue->rnext->Catom->Coordz,\
		 TempResidue->rnext->Oatom->Coordx,TempResidue->rnext->Oatom->Coordy,TempResidue->rnext->Oatom->Coordz,\
		 TempResidue->Catom->Coordx,TempResidue->Catom->Coordy,TempResidue->Catom->Coordz,\
		 TempResidue->Oatom->Coordx,TempResidue->Oatom->Coordy,TempResidue->Oatom->Coordz);
	}
}

void CalKappa()
{
	Residue *TempResidue;
	
	//前兩個residue賦值
	TestProtein.FirstResidue->Kappa = 360; 
	TestProtein.FirstResidue->rnext->Kappa = 360;
	
	//from first residue move to last residue
	for(TempResidue=TestProtein.FirstResidue;TempResidue->rnext->rnext->rnext->rnext;TempResidue=TempResidue->rnext)
	{
		//每一個胺基酸都要跟前後兩個比較(i-2,i,i+2)，但最後存取TempResidue->rnext->rnext
		//i = TempResidue->rnext->rnext ; i-2 = TempResidue ; i+2 = TempResidue->rnext->rnext->rnext->rnext ;
		TempResidue->rnext->rnext->Kappa = Angle\
		(TempResidue->rnext->rnext->CAatom->Coordx,TempResidue->rnext->rnext->CAatom->Coordy,TempResidue->rnext->rnext->CAatom->Coordz,\
		 TempResidue->CAatom->Coordx,TempResidue->CAatom->Coordy,TempResidue->CAatom->Coordz,\
		 TempResidue->rnext->rnext->rnext->rnext->CAatom->Coordx,TempResidue->rnext->rnext->rnext->rnext->CAatom->Coordy,TempResidue->rnext->rnext->rnext->rnext->CAatom->Coordz,\
		 TempResidue->rnext->rnext->CAatom->Coordx,TempResidue->rnext->rnext->CAatom->Coordy,TempResidue->rnext->rnext->CAatom->Coordz);	
	}
	
	//最後兩個residue賦值 
	TempResidue->rnext->rnext->Kappa = 360;
	TempResidue->rnext->rnext->rnext->Kappa = 360;
}

void CalAlpha()
{
	Residue *TempResidue;
	
	//前一個residue賦值
	TestProtein.FirstResidue->Alpha = 360;
	
	//from first residue move to last residue
	for(TempResidue=TestProtein.FirstResidue;TempResidue->rnext->rnext->rnext;TempResidue=TempResidue->rnext)
	{
		//每一個胺基酸都要跟前後兩個比較(i-1,i,i+1,i+2)，但最後存取TempResidue->rnext
		//i-1 = TempResidue ; i = TempResidue->rnext ; i+1 = TempResidue->rnext->rnext ; i+2 = TempResidue->rnext->rnext->rnext ;
		TempResidue->rnext->Alpha = TorsionAngle\
		(TempResidue->CAatom->Coordx,TempResidue->CAatom->Coordy,TempResidue->CAatom->Coordz,\
		 TempResidue->rnext->CAatom->Coordx,TempResidue->rnext->CAatom->Coordy,TempResidue->rnext->CAatom->Coordz,\
		 TempResidue->rnext->rnext->CAatom->Coordx,TempResidue->rnext->rnext->CAatom->Coordy,TempResidue->rnext->rnext->CAatom->Coordz,\
		 TempResidue->rnext->rnext->rnext->CAatom->Coordx,TempResidue->rnext->rnext->rnext->CAatom->Coordy,TempResidue->rnext->rnext->rnext->CAatom->Coordz);	
	}
	
	//最後兩個residue賦值 
	TempResidue->rnext->Alpha = 360;
	TempResidue->rnext->rnext->Alpha = 360;
}


void CalPhi()
{
	Residue *TempResidue;
	
	//前一個residue賦值
	TestProtein.FirstResidue->phi = 360;
	
	//from first residue move to last residue
	for(TempResidue=TestProtein.FirstResidue;TempResidue->rnext;TempResidue=TempResidue->rnext)
	{
		//每一個胺基酸都要跟前後兩個比較(i-1,i!!)，但最後存取TempResidue->rnext
		//i-1 = TempResidue ; i = TempResidue->rnext ;
		TempResidue->rnext->phi = TorsionAngle\
		(TempResidue->Catom->Coordx,TempResidue->Catom->Coordy,TempResidue->Catom->Coordz,\
		 TempResidue->rnext->Natom->Coordx,TempResidue->rnext->Natom->Coordy,TempResidue->rnext->Natom->Coordz,\
		 TempResidue->rnext->CAatom->Coordx,TempResidue->rnext->CAatom->Coordy,TempResidue->rnext->CAatom->Coordz,\
		 TempResidue->rnext->Catom->Coordx,TempResidue->rnext->Catom->Coordy,TempResidue->rnext->Catom->Coordz);	
	}

}

void CalPsi()
{
	Residue *TempResidue;
	
	//from first residue move to last residue
	for(TempResidue=TestProtein.FirstResidue;TempResidue->rnext;TempResidue=TempResidue->rnext)
	{
		//每一個胺基酸都要跟前後兩個比較(i,i+1!!)，但最後存取TempResidue
		//i = TempResidue ; i+1 = TempResidue->rnext ; 
		TempResidue->psi = TorsionAngle\
		(TempResidue->Natom->Coordx,TempResidue->Natom->Coordy,TempResidue->Natom->Coordz,\
		 TempResidue->CAatom->Coordx,TempResidue->CAatom->Coordy,TempResidue->CAatom->Coordz,\
		 TempResidue->Catom->Coordx,TempResidue->Catom->Coordy,TempResidue->Catom->Coordz,\
		 TempResidue->rnext->Natom->Coordx,TempResidue->rnext->Natom->Coordy,TempResidue->rnext->Natom->Coordz);	
	}
	
	//最後residue賦值 
	TempResidue->psi = 360;
}	
	
void OutputAngle()
{
	FILE *fptr2;
	fptr2 = fopen("./angle.txt","w");
	if(fptr2 == NULL)
	{
		printf ("Can't write to the file path: ./angle.txt \n"); 
		exit(0);
	}
	
	Residue *TempResidue;
	int i=0;
	fprintf(fptr2,"|Residue\t|TCO\t\t|Kappa\t\t|Alpha\t\t|Phi\t\t|Psi\n");
	for(TempResidue=TestProtein.FirstResidue;TempResidue->rnext;TempResidue=TempResidue->rnext)
	{
		fprintf(fptr2,"|%7d\t",i+1);
		fprintf(fptr2,"|%7.3f\t",TempResidue->TCO);
		fprintf(fptr2,"|%7.3f\t",TempResidue->Kappa);
		fprintf(fptr2,"|%7.3f\t",TempResidue->Alpha);
		fprintf(fptr2,"|%7.3f\t",TempResidue->phi);
		fprintf(fptr2,"|%7.3f\n",TempResidue->psi);
		i+=1;
	}	
	fclose(fptr2);
}
 

int main()
{
	//initialize
	TestProtein.FirstAtom = NULL;
	TestProtein.LastAtom = NULL;
	TestProtein.FirstResidue = NULL;
	TestProtein.LastResidue = NULL;
	
	ReadPDB("1crn");
	
//	printf("FirstResidueNum:%d\n",TestProtein.FirstResidue->Natom->ResNumber);
//	printf("LastResidueNum:%d\n",TestProtein.LastResidue->Natom->ResNumber);
	
	CalTCO();
	CalKappa();
	CalAlpha();
	CalPhi();
	CalPsi();

	OutputAngle(); 
//	DistanceMatrix(TestProtein.LastAtom->ResNumber);
	
	return 0;	
} 
