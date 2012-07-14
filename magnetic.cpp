#include <mpi.h>
#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <vector>
#include <time.h>
#include <fstream>
#include <algorithm>
#include <sstream>

using namespace std;

bool dumpData(string filename, int n, double *mn, double *mxyz, double *dm, double *hxyz)
{
	filebuf fb;
	fb.open (filename.c_str(), ios::out);
	ostream os(&fb);
	cout << "Dumping data to file " << filename << endl;
	for (int i = 0; i < n * 3; i += 3)
		os << mn[i] << " " << mn[i+1] << " " << mn[i+2] <<
		 " " << mxyz[i] << " " << mxyz[i+1] << " " << mxyz[i+2] <<
		 " " << hxyz[i] << " " << hxyz[i+1] << " " << hxyz[i+2] <<
		 " " << dm[i/3] << endl;
	fb.close();
}

void gener(float c,int sp,int sp1,int n, int N, double *mn, double *modl, double *cxyz,double *mxyz,double *vol,double *hc,double *hxyz,double *dm)
{
	vector <double> f;
	// elementy  f[0-d-1] funkcia raspredelinia polea vzaimodeistvia
	// elementy  f[d-d] funkcia raspredelinia polea vzaimodeistvia

	double kl=1E-4,kv=1E-12; //koefficienty dlay privedinia k santimentram
	double I0s=500.; //nachalnaya spontannay namagnichennost'

	//*****************************************************
	//creation of sample s zadannim raspedeleniem kriticheskih poley

	float sch;
	double pi=3.1415926,alfa,b, gamma;
	//zadadim koordinaty chastic i napravlyaushie cos
	srand(time(NULL)+7683347*time(NULL));
	for(int i=0; i<N; i+=3)
	{
		mn[i]=(double)rand()/(double)RAND_MAX;
		mn[i+1]=(double)rand()/(double)RAND_MAX;
		mn[i+2]=(double)rand()/(double)RAND_MAX;
		alfa=(double)rand()/(double)RAND_MAX*2*pi;
		b=(double)rand()/(double)RAND_MAX;
		gamma=(double)rand()/(double)RAND_MAX;
		if(gamma<0.5)b=-b;
		cxyz[i]=cos(alfa)*sqrt(1-b*b);
		cxyz[i+1]=sin(alfa)*sqrt(1-b*b);
		cxyz[i+2]=b;
	}
	float h0min = 20;
	float h0max = 55;
	//cout<< "\n Interval criticheskih poley h0min="<<h0min<<" h0max="<<h0max<<endl; 
	//raspredelim kriticheskie polya
	int ii=0;
	float d0;
	int	p= n/20;

	if (sp==1)
	{
		while(ii<p)
		{ 
			hc[ii]=20.+rand()%6; 
			ii++;
		}
		while(ii<p+2*p)
		{ 
			hc[ii]=25.+rand()%6; 
			ii++;
		}
		while(ii<p+2*p+4*p)
		{ 
			hc[ii]=30.+rand()%6; 
			ii++;
		}
		while(ii<p+2*p+4*p+6*p)
		{ 
			hc[ii]=35.+rand()%6; 
			ii++;
		}
		while(ii<p+2*p+4*p+6*p+4*p)
		{ 
			hc[ii]=40.+rand()%6; 
			ii++;
		}
		while(ii<p+2*p+4*p+6*p+4*p+2*p)
		{ 
			hc[ii]=45.+rand()%6; 
			ii++;
		}
		while(ii<p+2*p+4*p+6*p+4*p+2*p+p)
		{ 
			hc[ii]=50.+rand()%6; 
			ii++;
		}
	}

	//ravnomernoe raspredelenie kriticheskih poley
	if(sp==2)
	{
		for (int i=0; i<n; i++)
			hc[i]=(double)rand()/(double)RAND_MAX;	
	}

	// raspredelim diametry ispolzuem lognormalnoe raspredelenie
	ii=0;
	d0=0;
	if(sp1==3)
	{
		while(ii<p)
		{
			dm[ii]=0.3+(double)(rand()%1000)/10000.;
			ii++;
		}
		while(ii<p+5*p)
		{ 
			dm[ii]=0.4+(double)(rand()%1000)/10000.;
			ii++;
		}
		while(ii<p+5*p+5*p)//11
		{ 
			dm[ii]=0.5+(double)(rand()%1000)/10000.;  
			ii++;
		}
		while(ii<p+5*p+5*p+4*p)//15
		{ 
			dm[ii]=0.6+(double)(rand()%1000)/10000.; 
			ii++;
		}
		while(ii<p+5*p+5*p+4*p+3*p)//18
		{ 
			dm[ii]=0.7+(double)(rand()%1000)/10000.; 
			ii++;
		}
		while(ii<p+5*p+5*p+4*p+3*p+p)//19
		{ 
			dm[ii]=0.8+(double)(rand()%1000)/10000.; 
			ii++;
		}
		while(ii<p+5*p+5*p+4*p+3*p+p+p)//20
		{ 
			dm[ii]=0.9+(double)(rand()%1000)/10000.;
			ii++;
		}
	}

	//ravnomernoe raspredelenie, vse ravny 0.5
	if(sp1==2)
	{
		for (int i=0; i<n; i++)
			dm[i]=0.5;
	}
	ii=0;

	//Gauss
	if (sp1==1)
	{
		while(ii<p)
		{
			dm[ii]=0.3+(double)(rand()%1000)/10000.;
			ii++;
		}
		while(ii<p+2*p)
		{ 
			dm[ii]=0.4+(double)(rand()%1000)/10000.;
			ii++;
		}
		while(ii<p+2*p+4*p)//11
		{ 
			dm[ii]=0.5+(double)(rand()%1000)/10000.; 
			ii++;
		}
		while(ii<p+2*p+4*p+6*p)//15
		{ 
			dm[ii]=0.6+(double)(rand()%1000)/10000.;
			ii++;
		}
		while(ii<p+2*p+4*p+6*p+4*p)//18
		{ 
			dm[ii]=0.7+(double)(rand()%1000)/10000.;
			ii++;
		}
		while(ii<p+2*p+4*p+6*p+4*p+2*p)//19
		{ 
			dm[ii]=0.8+(double)(rand()%1000)/10000.;
			ii++;
		}
		while(ii<p+2*p+4*p+6*p+4*p+2*p+p)//20
		{ 
			dm[ii]=0.9+(double)(rand()%1000)/10000.; 
			ii++;
		}
	}

	float vob=0;
	for	(int i=0; i<n; i++)
	{
		vol[i]=(4./3.)*3.141*(dm[i]/2.)*(dm[i]/2.)*(dm[i]/2.);
		vob+= vol[i];
	}

	//ob'em i rasmery obrazca
	float	lObr = exp((1./3.)*log(vob/c));
	lObr = int(lObr + .5);

	//cout<<"\n lineiny rasmer obrazca="<<lObr<<"ob'em of sample "<< lObr*lObr*lObr<<endl;
	//cout<<"\n ob'em obrazca="<<vob<<endl;
	//cout<<"\n consentration="<<vob/lObr/lObr/lObr<<endl;

	//scorrectiruem coordinaty v sootvetstvie s lObr
	for(int i=0;i<N;i+=3)
	{
		mn[i]=mn[i]*lObr;
		mn[i+1]=mn[i+1]*lObr;
		mn[i+2]=mn[i+2]*lObr;
	}

	//dumpData("linear.txt", n, mn, mxyz, dm, hxyz);

	//zamenim peresekayuschiesya i blizkie chasticy
	int	a=0,  np=0, np1=0, it=0, itr=0;
    double  rmin=0, rx=0, ry=0, rz=0, r1=0; //chislo chastic do otbrasivania
    it=0;
    do
    {  
    	np=0;   
    	while (it<N)
    	{
    		for(int t=it+1; t<N; t+=3)
    		{
    			rx=fabs(mn[it]-mn[t]);
    			ry=fabs(mn[it+1]-mn[t+1]);
    			rz=fabs(mn[it+2]-mn[t+2]);
    			r1=sqrt((double)(rx*rx+ry*ry+rz*rz));
    			if (r1<(dm[it/3]+dm[t/3])/2. )
    			{
					//cout<<r1<<"   "<<dm[it/3]<<"   "<<dm[t/3]<<"   "<<(dm[it/3]+dm[t/3])/2.<<endl;
    				dm[it/3]=r1/2.; dm[t/3]=r1/2.;
					//cout<<r1<<"   "<<dm[it/3]<<"   "<<dm[t/3]<<"   "<<(dm[it/3]+dm[t/3])/2.<<endl;
    				np++;
	                //cout<<np<<endl;
    				np1++;
    			}
    		}
    		it+=3;
    		itr++;
    	}
    }while(np!=0);

    vob=0;
    for (int i=0; i<n; i++)
    {
    	vol[i]=(4./3.)*3.141*(dm[i]/2.)*(dm[i]/2.)*(dm[i]/2.);
    	vob+= vol[i];
    }

	//ob'em i rasmery obrazca
    lObr = exp((1./3.)*log(vob/c));
    lObr = int(lObr + .5);
	
	cout << "\nlinear size of sample = " << lObr << "\nvolume of sample " << lObr*lObr*lObr << endl;
	cout << "volume of particles = " << vob << endl;
	cout << "concentration = " << vob/lObr/lObr/lObr << endl;
	cout << "All " << np1 << " intersected particles separated" << endl;


	//podchitaem momenty
    for (int i=0; i<N; i+=3) 
    {
    	mxyz[i]= kv*vol[i/3]*cxyz[i]*I0s;
    	mxyz[i+1]= kv*vol[i/3]*cxyz[i+1]*I0s;
    	mxyz[i+2]= kv*vol[i/3]*cxyz[i+2]*I0s;
    }
	double hr=0, r2=0,r11=0, r111=0, mxr1=0, mxr2=0;//
	//pole formirovania obrazca

	double Hx=0.5, Hy=0., Hz=0, m0;
	//double Hx=0.5/sqrt(2.), Hy=0., Hz=-0.5/sqrt(2.), m0;
	cout<<"field of sidimentation Hsx="<<Hx<<"  Hsy="<<Hy<<"  Hsz="<<Hz<<endl;

	// set fields and scale to centimeters
	for( int i=0; i<N; i+=3)
	{
		hxyz[i] = Hx;
		hxyz[i+1] = Hy;
		hxyz[i+2] = Hz;

		mn[i] *= kl;
		mn[i+1] *= kl;
		mn[i+2] *= kl;

		dm[i/3] *= kl;
	}

	//dumpData("multiplied.txt", n, mn, mxyz, dm, hxyz);

	vector<double> ranjZ;
	vector<int> vspmZ;
	int nmaxf=-1, id=0, f1=0, prov=0;
	double maxf,maxf1;
	for(int i=0; i<N; i+=3)//Zapolnili massiv z sostovlyaushimi
	{
		ranjZ.push_back(mn[i+2]);
	}
	for(int i=0; i<n;i++)
		vspmZ.push_back(i);//zapolnili indexi ot 0 do n
	bool t=true;
	double buf1;
	int buf2,k=0;
	double hcmax=hc[0];
	while(t)
	{
		t=false;
		for(int j=0;j<ranjZ.size()-k-1;j++)
		{
			if(ranjZ[j]>ranjZ[j+1])
			{
				buf1=ranjZ[j];
				ranjZ[j]=ranjZ[j+1];
				ranjZ[j+1]=buf1;
				buf2=vspmZ[j];
				vspmZ[j]=vspmZ[j+1];
				vspmZ[j+1]=buf2;
				t=true;
			}
		}
		k++;
	}
	
	int i=3*vspmZ[0],j, flag=0;
	modl[i/3]=sqrt(hxyz[i]*hxyz[i]+hxyz[i+1]*hxyz[i+1]+hxyz[i+2]*hxyz[i+2]);
	m0= kv*vol[i/3]*I0s;
	cxyz[i]=hxyz[i]/modl[i/3];
	cxyz[i+1]=hxyz[i+1]/modl[i/3];
	cxyz[i+2]=hxyz[i+2]/modl[i/3];
	mxyz[i]=m0*cxyz[i];
	mxyz[i+1]=m0*cxyz[i+1];
	mxyz[i+2]=m0*cxyz[i+2];
	for(int j1=0; j1<n; j1++)
	{
		j=3*vspmZ[j1];
		if(i!=j)
		{
			r2= (mn[i]-mn[j])*(mn[i]-mn[j])+(mn[i+1]-mn[j+1])*(mn[i+1]-mn[j+1])+(mn[i+2]-mn[j+2])*(mn[i+2]-mn[j+2]);
			r11=sqrt(r2);
			mxr1= mxyz[i]*(mn[i]-mn[j])+mxyz[i+1]*(mn[i+1]-mn[j+1])+mxyz[i+2]*(mn[i+2]-mn[j+2]);
			r111=(r11*r11*r11*r11*r11);
			r11=r11*r11*r11;
			hxyz[j]+=3.*mxr1*(mn[i]-mn[j])/r111-mxyz[i]/r11;
			hxyz[j+1]+=3.*mxr1*(mn[i+1]-mn[j+1])/r111-mxyz[i+1]/r11;
			hxyz[j+2]+=3.*mxr1*(mn[i+2]-mn[j+2])/r111-mxyz[i+2]/r11;
		}
	}
	for (int i1=1; i1<n; i1++)
	{
		//if(i1%100==0)cout<<i1<<endl;
		i=3*vspmZ[i1];
		modl[i/3]=sqrt(hxyz[i]*hxyz[i]+hxyz[i+1]*hxyz[i+1]+hxyz[i+2]*hxyz[i+2]);
		m0=kv*vol[i/3]*I0s;
		cxyz[i]=hxyz[i]/modl[i/3];
		cxyz[i+1]=hxyz[i+1]/modl[i/3];
		cxyz[i+2]=hxyz[i+2]/modl[i/3];
		mxyz[i]=m0*cxyz[i];	
		mxyz[i+1]=m0*cxyz[i+1];
		mxyz[i+2]=m0*cxyz[i+2];
		for (int j1=0; j1<n; j1++)//n
		{
			j=3*vspmZ[j1];
			if(i!=j)
			{
				r2=(mn[i]-mn[j])*(mn[i]-mn[j])+(mn[i+1]-mn[j+1])*(mn[i+1]-mn[j+1])+(mn[i+2]-mn[j+2])*(mn[i+2]-mn[j+2]);
				r11=sqrt(r2);
				mxr1= mxyz[i]*(mn[i]-mn[j])+mxyz[i+1]*(mn[i+1]-mn[j+1])+mxyz[i+2]*(mn[i+2]-mn[j+2]);
				r111=(r11*r11*r11*r11*r11);
				r11=r11*r11*r11;      
				hxyz[j]+=3.*mxr1*(mn[i]-mn[j])/r111-mxyz[i]/r11;
				hxyz[j+1]+=3.*mxr1*(mn[i+1]-mn[j+1])/r111-mxyz[i+1]/r11;
				hxyz[j+2]+=3.*mxr1*(mn[i+2]-mn[j+2])/r111-mxyz[i+2]/r11;
			}
		}
		flag=1;
		do
		{   
			flag=0;maxf1=0;
			for(int it=0; it<=i1; it++)
			{
				id=3*vspmZ[it];
				modl[id/3]=sqrt(hxyz[id]*hxyz[id]+hxyz[id+1]*hxyz[id+1]+hxyz[id+2]*hxyz[id+2]);
				if(hcmax<hc[id/3])
				{
					hcmax=hc[id/3];                 
				}
				if(mxyz[id]*hxyz[id]+mxyz[id+1]*hxyz[id+1]+mxyz[id+2]*hxyz[id+2]<0 && modl[id/3]>hc[id/3])
				{
					flag++;
					maxf=modl[id/3];
					if(maxf1<=maxf)
					{
						nmaxf=id; maxf1=maxf;
					}   
				}
			}
			if(nmaxf!=-1 && flag!=0 && mxyz[nmaxf]*hxyz[nmaxf]+mxyz[nmaxf+1]*hxyz[nmaxf+1]+mxyz[nmaxf+2]*hxyz[nmaxf+2]<0)
			{
				for (int j1=0; j1<n; j1++)//n
				{
					j=3*vspmZ[j1];
					if(nmaxf!=j)
					{
						//raschitaem resultiruyuchee pole sozdannoe chasticei 0 na ostalnih chacticah
						r2=(mn[nmaxf]-mn[j])*(mn[nmaxf]-mn[j]) + (mn[nmaxf+1]-mn[j+1])* (mn[nmaxf+1]-mn[j+1])+(mn[nmaxf+2]-mn[j+2])* (mn[nmaxf+2]-mn[j+2]);
						r11=sqrt(r2);
						mxr1= mxyz[nmaxf]*(mn[nmaxf]-mn[j])+mxyz[nmaxf+1]*(mn[nmaxf+1]-mn[j+1])+mxyz[nmaxf+2]*(mn[nmaxf+2]-mn[j+2]);
						r111=(r11*r11*r11*r11*r11);
						r11=r11*r11*r11;
						hxyz[j]-=3.*mxr1*(mn[nmaxf]-mn[j])/r111-mxyz[nmaxf]/r11;
						hxyz[j+1]-=3.*mxr1*(mn[nmaxf+1]-mn[j+1])/r111-mxyz[nmaxf+1]/r11;
						hxyz[j+2]-=3.*mxr1*(mn[nmaxf+2]-mn[j+2])/r111-mxyz[nmaxf+2]/r11;
						mxyz[nmaxf]=-mxyz[nmaxf];
						mxyz[nmaxf+1]=-mxyz[nmaxf+1];
						mxyz[nmaxf+2]=-mxyz[nmaxf+2];
						mxr1= mxyz[nmaxf]*(mn[nmaxf]-mn[j])+mxyz[nmaxf+1]*(mn[nmaxf+1]-mn[j+1])+mxyz[nmaxf+2]*(mn[nmaxf+2]-mn[j+2]);
						hxyz[j]+=3.*mxr1*(mn[nmaxf]-mn[j])/r111-mxyz[nmaxf]/r11;
						hxyz[j+1]+=3.*mxr1*(mn[nmaxf+1]-mn[j+1])/r111-mxyz[nmaxf+1]/r11;
						hxyz[j+2]+=3.*mxr1*(mn[nmaxf+2]-mn[j+2])/r111-mxyz[nmaxf+2]/r11;
					}
				}
			}
		}while(flag!=0);
	}
}

int main (int argc, char* argv[])
{
	int rank, size;
	MPI_Status Status;
	MPI_Request req;
	MPI_Init (&argc, &argv);
	MPI_Comm_rank (MPI_COMM_WORLD, &rank);
	MPI_Comm_size (MPI_COMM_WORLD, &size);

	int n = 10000;				//Kolichestvo chastic
	if (rank==0)
		cout << "number of particles: n = " << n << endl;
	int state=1;				//Initial state 1-NS, 2 - Orient, 3- ori X, 4-ori2 z, 5 - ori2
	int sp=1;					//Spectr criticheskyh poley  1- Gauss, 2- Ravnomerny
	int sp1=3;					//Raspredelenie ob'emov 1- Gauss, 2- Ravnomerny, 3 - Lognormalny
	float hstp=1;				//Step of field
	float c=0.05/100.;			//Consentration
	int N=3*n;
	double Hx, Hy, Hz;    
	double *mn = new double [N];//koordinaty chastic x,y,z
	double *modl=new double [n];//module  polya na chastice
	double *cxyz=new double [N];//cosx,cosy,cosz
	double *mxyz=new double [N];//magnitnie momenty, mx,my,mz
	double *vol=new double [n]; //volume
	double *hc=new double [n];	//critical field
	double *hxyz=new double [N];//summarnie komponenty polya na zadannoi chastice
	double *dm=new double [n];	//diameter
	vector<int> nr;				//nomera neravnovvesnih chastic na tekushem shage
	int *per;					//chastici dlya perevorota
	int kper;					//kolichestvo chastic dlya perevorota
	int nnc=1;					//kolvo neravnovesnih
	int shag=0;					//nomer shaga razmagnichivaniya
	if (rank==0)
		gener(c, sp, sp1, n, N, mn, modl, cxyz, mxyz, vol, hc, hxyz, dm);//Generiruem dannie
	MPI_Bcast(mn,N,MPI_DOUBLE,0,MPI_COMM_WORLD);//Rassilaem na vseh
	MPI_Bcast(modl,n,MPI_DOUBLE,0,MPI_COMM_WORLD);//Rassilaem na vseh
	MPI_Bcast(cxyz,N,MPI_DOUBLE,0,MPI_COMM_WORLD);//Rassilaem na vseh
	MPI_Bcast(mxyz,N,MPI_DOUBLE,0,MPI_COMM_WORLD);//Rassilaem na vseh
	MPI_Bcast(vol,n,MPI_DOUBLE,0,MPI_COMM_WORLD);//Rassilaem na vseh
	MPI_Bcast(hc,n,MPI_DOUBLE,0,MPI_COMM_WORLD);//Rassilaem na vseh
	MPI_Bcast(hxyz,N,MPI_DOUBLE,0,MPI_COMM_WORLD);//Rassilaem na vseh
	MPI_Bcast(dm,n,MPI_DOUBLE,0,MPI_COMM_WORLD);//Rassilaem na vseh
	double I0s=500.;
	float	IrNsZ=0, IrNsX=0, IrNsY=0, IrsX=0, IrsY=0, IrsZ=0;
	for(int i=0; i<N; i+=3)
	{
		IrNsZ = IrNsZ + I0s*vol[i/3]* cxyz[i+2];
		IrNsY = IrNsY + I0s*vol[i/3]* cxyz[i+1];
		IrNsX = IrNsX + I0s*vol[i/3]* cxyz[i];
		IrsZ = IrsZ + I0s*vol[i/3]* fabs(cxyz[i+2]); 
		IrsY = IrsY + I0s*vol[i/3]* fabs(cxyz[i+1]);
		IrsX = IrsX + I0s*vol[i/3]* fabs(cxyz[i]);
	}
	if (rank==0)
		cout<<" Initial state of sample"<<endl;
	if (rank==0)
		cout<<"IrNsX="<<IrNsX<<"\tIrNsY="<<IrNsY<<"\tIrNsZ="<<IrNsZ<<endl;
	if (rank==0)
		cout<<"IrsX="<<IrsX<<"\tIrsY="<<IrsY<<"\tIrsZ="<<IrsZ<<endl;

	// dump mn and mxyz to file
	if (rank == 0)
		dumpData("initial.txt", n, mn, mxyz, dm, hxyz);

	int step, znak=-1, zz=0;//peremennie neobhodimie dlya izmeneniya polya
	double vich;
	do{
		//naidem samoe bolshoe pole
		if(rank==0)cout<<"******************** Itteration "<<zz+1<<endl;
		if(zz==0)
		{
			float hx=300, hy=300, hz=300;
			Hx=hxyz[0];
			Hy=hxyz[1];
			Hz=hxyz[2];
			for(int i=0; i<N; i+=3)
			{
				if((Hx*Hx+Hy*Hy+Hz*Hz)<(hxyz[i]*hxyz[i]+hxyz[i+1]*hxyz[i+1]+hxyz[i+2]*hxyz[i+2]))
				{
					Hx=hxyz[i];
					Hy=hxyz[i+1];
					Hz=hxyz[i+2];
				}
			}
			Hx=fabs(Hx);
			Hy=fabs(Hy);
			Hz=fabs(Hz);
			Hx+=hx;
			Hy+=hy;
			Hz+=hz;
			Hx=sqrt(Hx*Hx+Hy*Hy+Hz*Hz);
			Hy=0;
			Hz=0;
		}
		for(int i=0; i<N; i+=3)
		{
			hxyz[i]+=Hx;
			hxyz[i+1]+=Hy;
			hxyz[i+2]+=Hz;
			modl[i/3]=sqrt(hxyz[i]*hxyz[i]+hxyz[i+1]*hxyz[i+1]+hxyz[i+2]*hxyz[i+2]);   
		}
		if(rank==0)
			cout<<"maximal module of internal field Hx="<<Hx<<"  Hy="<<Hy<<"  Hz="<<Hz<<endl;
		do
		{
			nnc=0;
			int id;
			nr.erase(nr.begin(),nr.end());
			for(int i=0;i<n;i++)//naidem neravnovesnie i zapishem nomera
			{
				id=3*i;
				if(mxyz[id]*hxyz[id]+mxyz[id+1]*hxyz[id+1]+mxyz[id+2]*hxyz[id+2]<0 && modl[id/3]> hc[id/3])
				{
					nnc++;
					nr.push_back(id);
				}   
			}
			if(nnc!=0)
			{
				for(int i=0;i<nnc;i++)
				{
					id=nr[i];
					mxyz[id]*=(-1);
					mxyz[id+1]*=(-1);
					mxyz[id+2]*=(-1);
					cxyz[id]*=(-1);
					cxyz[id+1]*=(-1);
					cxyz[id+2]*=(-1);
				}
				//Raschitaem energiu
				float E=0, r;
				int id1,id2;
				int x=n/size;
				int low=rank*x;
				int high=low+x;
				if(rank==size-1)high=n;
				for(int i=low;i<high;i++)
				{
					for(int j=i+1;j<n;j++)
					{
						id1=3*i, id2=3*j;
						r=sqrt((mn[id1]-mn[id2])*(mn[id1]-mn[id2])+(mn[id1+1]-mn[id2+1])*(mn[id1+1]-mn[id2+1])+(mn[id1+2]-mn[id2+2])*(mn[id1+2]-mn[id2+2]));
						E-=((mxyz[id1]*mxyz[id2]+mxyz[id1+1]*mxyz[id2+1]+mxyz[id1+2]*mxyz[id2+2])/(r*r*r)-3*(mxyz[id1]*(mn[id1]-mn[id2])+mxyz[id1+1]*(mn[id1+1]-mn[id2+1])+mxyz[id1+2]*(mn[id1+2]-mn[id2+2]))*(mxyz[id2]*(mn[id2]-mn[id1])+mxyz[id2+1]*(mn[id2+1]-mn[id1+1])+mxyz[id2+2]*(mn[id2+2]-mn[id1+2]))/(r*r*r*r*r));
					}
				}
				float E_r;
					//cout<<E<<endl;
				MPI_Reduce(&E, &E_r, 1, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
				E=E_r;
				if(rank==0)cout<<"E="<<E<<"\t"<<"  number of nonequlibrium="<<nnc<<endl;
				float r2,r11,r111,mxr1,mxr2;
				x=n/size;
				low=rank*x;
				high=low+x;
				if(rank==size-1)high=n;
				for(int j=low;j<high;j++)   
				{ 
					id2=3*j;
					hxyz[id2]=Hx;
					hxyz[id2+1]=Hy;
					hxyz[id2+2]=Hz;
					for(int i=j+1;i<n;i++)
					{
						id1=3*i;
						r2= (mn[id1]-mn[id2])*(mn[id1]-mn[id2])+(mn[id1+1]-mn[id2+1])*(mn[id1+1]-mn[id2+1])+(mn[id1+2]-mn[id2+2])*(mn[id1+2]-mn[id2+2]);
						r11=sqrt(r2);
						r111=r11*r11*r11*r11*r11;
						r11=r11*r11*r11;
						mxr1= mxyz[id1]*(mn[id1]-mn[id2])+mxyz[id1+1]*(mn[id1+1]-mn[id2+1])+mxyz[id1+2]*(mn[id1+2]-mn[id2+2]);
						hxyz[id2]+=3.*mxr1*(mn[id1]-mn[id2])/r111-mxyz[id1]/r11;
						hxyz[id2+1]+=3.*mxr1*(mn[id1+1]-mn[id2+1])/r111-mxyz[id1+1]/r11;
						hxyz[id2+2]+=3.*mxr1*(mn[id1+2]-mn[id2+2])/r111-mxyz[id1+2]/r11;
						mxr2= mxyz[id2]*(mn[id2]-mn[id1])+mxyz[id2+1]*(mn[id2+1]-mn[id1+1])+mxyz[id2+2]*(mn[id2+2]-mn[id1+2]);
						hxyz[id1]+=3.*mxr2*(mn[id2]-mn[id1])/r111-mxyz[id2]/r11;
						hxyz[id1+1]+=3.*mxr2*(mn[id2+1]-mn[id1+1])/r111-mxyz[id2+1]/r11;
						hxyz[id1+2]+=3.*mxr2*(mn[id2+2]-mn[id1+2])/r111-mxyz[id2+2]/r11;
					}
					modl[j]=sqrt(hxyz[j]*hxyz[j]+hxyz[j+1]*hxyz[j+1]+hxyz[j+2]*hxyz[j+2]);
				}
				for(int i=0;i<n;i++)
				{
					if((i<low)||(i>=high))
					{
						id2=i*3;
						hxyz[id2]=0;
						hxyz[id2+1]=0;
						hxyz[id2+2]=0;
						modl[i]=0;
					}
				}
				double *hxyz_r=new double [N];
				double *modl_r=new double [n];
				MPI_Reduce(hxyz, hxyz_r, N, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
				MPI_Reduce(modl, modl_r, n, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
				for(int i=0;i<n;i++)
				{
					id2=3*i;
					hxyz[id2]=hxyz_r[id2];
					hxyz[id2+1]=hxyz_r[id2+1];
					hxyz[id2+2]=hxyz_r[id2+2];
					modl[i]=modl_r[i];
				}
				delete [] hxyz_r;
				delete [] modl_r;
				MPI_Bcast(modl,n,MPI_DOUBLE,0,MPI_COMM_WORLD);
				MPI_Bcast(hxyz,N,MPI_DOUBLE,0,MPI_COMM_WORLD);
				double I0s=500.;
				float	IrNsZ=0, IrNsX=0, IrNsY=0, IrsX=0, IrsY=0, IrsZ=0, Mx=0, My=0, Mz=0;
				for(int i=0; i<N; i+=3)
				{
					IrNsZ = IrNsZ + I0s*vol[i/3]* cxyz[i+2];
					IrNsY = IrNsY + I0s*vol[i/3]* cxyz[i+1];
					IrNsX = IrNsX + I0s*vol[i/3]* cxyz[i];
					IrsZ = IrsZ + I0s*vol[i/3]* fabs(cxyz[i+2]); 
					IrsY = IrsY + I0s*vol[i/3]* fabs(cxyz[i+1]);
					IrsX = IrsX + I0s*vol[i/3]* fabs(cxyz[i]);
					Mx+=mxyz[i];
					My+=mxyz[i+1];
					Mz+=mxyz[i+2];
				}
				if(rank==0) {
					cout<<"IrNsX="<<IrNsX<<"\tIrNsY="<<IrNsY<<"\tIrNsZ="<<IrNsZ<<endl;

					ostringstream filename;
					filename << "after-" << zz << ".txt";
					//dumpData(filename.str(), n, mn, mxyz, dm, hxyz);
				}
			}
		}while(nnc!=0);
		Hx=fabs(Hx);
		if(Hx>10)
		{
			ostringstream ost;
			ost << (int)Hx;
			step=ost.str().size()-1;
			vich=exp(step * log(10.));
			if((Hx-vich)<(vich))vich=vich/10;
		}
		else
		{
			vich=0.1;
		}
		Hx-=vich;
		Hx*=znak;
		znak*=-1;
		zz++;
	}while(fabs(Hx)>0.1);
	if (rank == 0)
		dumpData("final.txt", n, mn, mxyz, dm, hxyz);
	MPI_Barrier(MPI_COMM_WORLD);//Sinhroniziruemsya
	MPI_Finalize();
	return 0;
}
