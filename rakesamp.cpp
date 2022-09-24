#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include <vector>
#include <stdlib.h>
#include <time.h>
#include <stdio.h>
#include <string.h>
#include <cstring>
#include <list>
#include <sstream>



// This is an implementation of the RakeSamp algorithm for RNA secondary structures without pseudoknots;
// prepared by Ricky X. F. Chen (School of Mathematics, Hefei University of Technology), chenshu731@sina.com
// webpage: http://maths.hfut.edu.cn/2020/1202/c7700a242703/page.htm

using namespace std;

vector<int> git(list<vector<int> >&,int);
int maxguy(vector<int>);
vector<int> generateRandom(vector<int> &);
vector<int> firstk(vector<int>, int);

int finesampz_par(vector<int>&,vector<int>&, vector<int>&, vector<int>&, vector<int>&, list<vector<int>>&, string);
vector<int> rake2dot(vector<int> , vector<int> , int , int , int );
vector<int> pairroot(list<vector<int>> , list<vector<int>> , int , int , int , int , int);
vector<int> con_nest(list<vector<int>> , vector<int> , int , int , int , int);
void outdot(vector<int> , string );



//this block is the function that returns the j-th vector in the list of vectors, j=1 gives the very first vector;
vector<int> git(list<vector<int>>&l, int j){
    list<vector<int>>::iterator it = l.begin();
    for(int i=0; i<j-1; i++){
        ++it;
    }
    return *it;
}
//end of block




//this block is the function that computes the maximum element in the passed vector
int maxguy(vector<int> A){
	int n = *max_element(A.begin(),A.end());
	return n;
}
//end of block




 
 
// Function to generate a random permutation of the passed sequence
vector<int> generateRandom(vector<int> &v)
{
   // vector<int> v(n);
    vector<int> vout;
 
 while(v.size()){
 
   // Size of the vector
    int n = v.size();
 
    // Generate a random number
    srand(time(NULL));
 
    // Make sure the number is within
    // the index range
    int index = rand() % n;
 
    // Get random number from the vector
    int num = v[index];
 
    // Remove the number from the vector
    swap(v[index], v[n - 1]);
    v.pop_back();
    

    
        vout.push_back(num);
    
}

    return vout;
}

// function returns the subsequence containing the first k elements
vector<int> firstk(vector<int> v, int k){
	int n=v.size();
	vector<int> kout;
	for(int i=0;i<k;i++){
		kout.push_back(v[i]);
	}
	return kout;
}


// getting paramters from a txt file
int finesampz_par(vector<int>&ex,vector<int>&hz, vector<int>&p, vector<int>&g, vector<int>&t, list<vector<int>>&m_ij, string netname){
	// in the file netname.txt, (1) the first row specifies exterior loop ed, el, separated by a space; 
	// (2) the 2nd row the number of helices of size i>=1; 
	//(3) the 3rd row, the number of hairpins of length i>=0; 
	// (4) the 4th row, the number of bulges of length i>=0; (5) the 5th row, the number of
	// interior loops of length i>= 0, (6) starting from the 6th row, the number of multiloops of length i>=0, where the first entry in the row
	// is the degree of the multiloop;
 // the parameters in the file should be compatible, and the code will not check the compatibility
 
 
	 ex.resize(0);
	 hz.resize(0);
	 p.resize(0);
	 g.resize(0);
	 t.resize(0);
	 
	 vector<int> triplet;
	 int type=1; 
	 
  string line3;
  ifstream my3file (netname+".txt");
  if (my3file.is_open())
  {
    while ( getline (my3file,line3) )
    {
     
 	char * cstr = new char [line3.length()+1];
  std::strcpy (cstr, line3.c_str());

  // cstr now contains a c-string copy of line3

	char * pp = std::strtok (cstr," ");
  triplet.push_back(atoi(pp));
	while (pp!=0)
  {
    pp = std::strtok(NULL," ");
    triplet.push_back(atoi(pp));
    
  }
  
  // the numbers in a row from the file will be stored in the vector triplet;
		  if(type==1){
		  	for(int i=0; i< triplet.size(); i++){
		  		ex.push_back(triplet[i]);
			  }
			  type++;
		  }
		  
		  else if(type==2){
		  	for(int i=0; i< triplet.size(); i++){
		  		hz.push_back(triplet[i]);
			  }
			  type++;
		  }
		  
		  else if(type==3){
		  	for(int i=0; i< triplet.size(); i++){
		  		p.push_back(triplet[i]);
			  }
			  type++;
		  }
		  
		  else if(type==4){
		  	for(int i=0; i< triplet.size(); i++){
		  		g.push_back(triplet[i]);
			  }
			  type++;
		  }
		  
		 else  if(type==5){
		  	for(int i=0; i< triplet.size(); i++){
		  		t.push_back(triplet[i]);
			  }
			  type++;
		  }
		  
		  
		  else if(type>5){
		  	m_ij.push_back(triplet);
			  type++;
		  }
  
	
	triplet.resize(0);

  delete[] cstr;
 
    }
    my3file.close();
  }
  else cout << "Unable to open file"; 
//end of block

return (type-6);  //return the number of multiloops specified in the file
}
//end of block


// rake2dot translate a rake into a sequence of -1, 1, 0, and numbers larger than b+k+1
vector<int> rake2dot(vector<int> S, vector<int> L, int b, int k, int h){
	vector<int> dotseq;
	for(int i=0; i< S.size(); i++){
		dotseq.push_back(-1);
	}
	for(int i=0; i< L.size(); i++){
		if(L[i]<b+k+2){
			dotseq.push_back(0);
		}
		else {
			dotseq.push_back(L[i]);
		}
	}
	for(int i=0; i<S.size(); i++){
		dotseq.push_back(1);
	}
	
	return dotseq;
}


// Note that each marked label, i.e., those larger than b+k+1, is associated to exact one unmarked label, i.e., those smaller than
// b+k+1, that is the root of a rake. The function returns the indices in the sequences of rakes corresponding to those found unmarked label.
vector<int> pairroot(list<vector<int>> S, list<vector<int>> L, int b, int k, int h, int ed, int el){
	
	vector<int>  rootindex;   //recorded starting from b+k+2; position starting from 1 due to git function;
//rootindex.resize(0);
vector<int>  temp, tempS;

vector<int> used;
int minroot, minindex; 


if((ed>2)||((ed==2)&&(el>0))){

minroot=b+k+1;
for(int j=1; j<h+1; j++){
	temp.resize(0);
	temp=git(L,j);
	if(maxguy(temp)<(b+k+2)){
			tempS.resize(0);
			tempS=git(S,j);
			if(tempS[0]<minroot){
				minroot=tempS[0];
				minindex=j;
			}
			
			
		
			}
		}

used.push_back(minroot);
rootindex.push_back(minindex);

for(int i=1; i<h; i++){
	minroot=b+k+1;
	for(int j=1; j<h+1; j++){
		temp.resize(0);
		temp=git(L,j);
		
		if(maxguy(temp)<(b+k+i+2)){
			tempS.resize(0);
			tempS=git(S,j);
			for(int q=0; q<used.size(); q++){
				if(tempS[0]==used[q]){
					break;
				}
				else if((q==(used.size()-1))&&(tempS[0]<minroot)){
					minroot=tempS[0];
					minindex=j;
				}
			}
		}
	}
	used.push_back(minroot);
	rootindex.push_back(minindex);
}
}

else{
minroot=b+k+1;
for(int j=1; j<h+1; j++){
	temp.resize(0);
	temp=git(L,j);
	if(maxguy(temp)<(b+k+2)){
			tempS.resize(0);
			tempS=git(S,j);
			if(tempS[0]<minroot){
				minroot=tempS[0];
				minindex=j;
			}
			
			
		
			}
		}

used.push_back(minroot);
rootindex.push_back(minindex);

for(int i=1; i<h-1; i++){ //there is one rake less for this case
	minroot=b+k+1;
	for(int j=1; j<h+1; j++){
		temp.resize(0);
		temp=git(L,j);
		
		if(maxguy(temp)<(b+k+i+2)){
			tempS.resize(0);
			tempS=git(S,j);
			for(int q=0; q<used.size(); q++){
				if(tempS[0]==used[q]){
					break;
				}
				else if((q==(used.size()-1))&&(tempS[0]<minroot)){
					minroot=tempS[0];
					minindex=j;
				}
			}
		}
	}
	used.push_back(minroot);
	rootindex.push_back(minindex);
}	
}


return rootindex;

} 



// concatenation and nesting of the -1-0-1-marked expression of the rakes, where each marked label
// is replaced by the expression of the rake with the index corresponding to the marked label;
// evetually, all marked labels will be eliminated.
vector<int> con_nest(list<vector<int>> Sdot, vector<int> pair, int b, int k, int h, int posroot){
			
vector<int> temp1, temp2, outseq;
int x;
temp1.resize(0);
temp2.resize(0);

temp2=git(Sdot, posroot); // posroot is the index of the rake where b+k+1 locates


while(maxguy(temp2)>1){
	outseq.resize(0);
	for(int i=0; i< temp2.size(); i++){
		if(temp2[i]<2){
			outseq.push_back(temp2[i]);
		}
		else{
			temp1.resize(0);
			x=temp2[i]-b-k-1;
			temp1=git(Sdot, pair[x-1]);
			for(int j=0; j< temp1.size(); j++){
				outseq.push_back(temp1[j]);
			}
		}
	}
	
	temp2.resize(0);
	for(int i=0; i<outseq.size(); i++){
		temp2.push_back(outseq[i]);
	}
}

return outseq;


}


// transform the -1-0-1 sequence into bracket-dot expression written into a txt file.
void outdot(vector<int> connest, string netname){
	ofstream stru;
	stru.open(netname+"_structure.txt");
	
	for(int i=1; i< connest.size()-1; i++){ // remove the auxilioury big arc
	if(connest[i]==-1){
		stru<<"(";
	}
	
	else if(connest[i]==1){
		stru<<")";
	}
	
	else if(connest[i]==0){
		stru<<".";
	}
	
}

stru.close();
}



//this block is the main block
int main (int argc, char** argv){
	

string netname; // netname for parameters and output secondary structures
//netname="exam";
cin>> netname;
//cin>> netname1;
//netname1="test3";

int b=0, k=0, h=0, ed, el;
//int P, G, T, M;
vector<int> p, g, t, hz, ex;
list< vector<int> > m_ij;

int nmloop=finesampz_par(ex, hz, p, g, t, m_ij, netname); // number of multiloops

ed=ex[0], el=ex[1];
for(int i=0; i<hz.size(); i++){
	h=h+hz[i];
}

for(int i=0; i<hz.size(); i++){
	b=b+(i+1)*hz[i];
}

k=k+el;
for(int i=0; i<p.size(); i++){
	k=k+i*p[i];
}

for(int i=0; i<g.size(); i++){
	k=k+i*g[i];
}

for(int i=0; i<t.size(); i++){
	k=k+i*t[i];
}

vector<int> m_j;
for(int i=1; i<nmloop+1; i++){
	m_j.resize(0);
	m_j=git(m_ij,i);
	for(int j=1; j<m_j.size(); j++){
		k=k+(j-1)*m_j[j];
	}
}

// end of collecting paramters


if((ed==2)&&(el==0)&&(h==1)){  //the trivial case h=1;
	
		ofstream stru;
	stru.open(netname+"_structure.txt");
	
	for(int i=0; i< b; i++){
	
		stru<<"(";
	}
	for(int i=0; i< k; i++){
	
		stru<<".";
	}
	for(int i=0; i< b; i++){
	
		stru<<")";
	}


stru.close();

}

else{

vector<int> X,Y;

for(int i=1;i<b+k+1;i++){
	Y.push_back(i);
}

for(int i=0;i<h-1;i++){
	X.push_back(b+k+2+i);
}


vector<int> Lev, Mlv, YY, temp;

int posroot, counter=0; // counter counting the number of rakes push back into the list

YY=generateRandom(Y);
Lev=firstk(YY,k);
Mlv=generateRandom(X);

list<vector<int> > S, L;
int posl=0, posm=0; // position in Lev and Mlv, resp.
int rnd;

for(int i=0;i<p.size();i++){
	if(p[i]!=0){
		for(int j=0;j<p[i];j++){
			temp.resize(0);
			for(int pp=0;pp<i;pp++){
				temp.push_back(Lev[posl]);
				posl++;
			}
			L.push_back(temp);
			counter++;
		}
	}
}
// completion of hairpins

srand(time(NULL));

for(int i=0;i<g.size();i++){
	if(g[i]!=0){
		for(int j=0;j<g[i];j++){
			temp.resize(0);
			rnd=rand() % 2;
			if(rnd==0){
				temp.push_back(Mlv[posm]);
				if(Mlv[posm]==(b+k+h)){
					posroot=counter+1;
				}
				posm++;
				for(int pp=0;pp<i;pp++){
				temp.push_back(Lev[posl]);
				posl++;
				}
			}
			else{
				for(int pp=0;pp<i;pp++){
				temp.push_back(Lev[posl]);
				posl++;
				}
				temp.push_back(Mlv[posm]);
				if(Mlv[posm]==(b+k+h)){
					posroot=counter+1;
				}
				posm++;
			}
			
			L.push_back(temp);
			counter++;
		}
	}
}
// completion of bulges

srand(time(NULL));

for(int i=0;i<t.size();i++){
	if(t[i]!=0){
		for(int j=0;j<t[i];j++){
			temp.resize(0);
			rnd=rand() % (i-1);
			for(int pp=0;pp<rnd+1;pp++){
				temp.push_back(Lev[posl]);
				posl++;
				}
				temp.push_back(Mlv[posm]);
				if(Mlv[posm]==(b+k+h)){
					posroot=counter+1;
				}
				posm++;
				for(int pp=rnd+2;pp<i+1;pp++){
				temp.push_back(Lev[posl]);
				posl++;
				}
			
			
			L.push_back(temp);
			counter++;
		}
	}
}

//completion of interior loops;

int mult=1;
vector<int> multls, mulpos, msize, Emv;

//multls=git(m_ij,mult);

while(mult<nmloop+1){
	multls=git(m_ij, mult);
	for(int i=1;i<multls.size();i++){
		
	if(multls[i]!=0){
		
	
		for(int j=0;j<multls[i];j++){
			msize.resize(0);
			temp.resize(0);
			for(int qq=0;qq<(i-1+multls[0]-1);qq++){
			msize.push_back(qq); // note: the sequence starting from 0;
			}
			mulpos=firstk(generateRandom(msize),multls[0]-1);
			
			for(int pp=0;pp<(i-1+multls[0]-1);pp++){
				for(int tt=0;tt<(multls[0]-1);tt++){
				
				if((pp==mulpos[tt])){
					temp.push_back(Mlv[posm]);
					if(Mlv[posm]==(b+k+h)){
					posroot=counter+1;
				}
					posm++;
					break;
				}
				else if((pp!=mulpos[tt])&&(tt==(multls[0]-2))){
					temp.push_back(Lev[posl]);
					posl++;
				}
			}
				
				}
		
			
			
			L.push_back(temp);
			counter++;
			mulpos.resize(0);
		}
	}
	
}

mult++;
multls.resize(0);
//multls=git(m_ij,mult);
}
//completion of multiloops
//rnd=rand() % (ed-1);

if((ed>2)||((ed==2)&&(el>0))){

posroot=h+1;

if(ed==2){
	temp.resize(0);
	for(int kk=0; kk<el; kk++){
		temp.push_back(Lev[posl]);
		posl++;
	}
}

if(ed>2){

mulpos.resize(0);
//for(j=0;j<multls[i];j++){
			msize.resize(0);
			temp.resize(0);
			for(int qq=0;qq<(ed-2+el);qq++){
			msize.push_back(qq);
			}

			mulpos=firstk(generateRandom(msize),ed-2); // ed-2>0 required here
			
			for(int pp=0;pp<(ed-2+el);pp++){
				for(int tt=0;tt<(ed-2);tt++){  // ed-2>0 required here
				
				if((pp==mulpos[tt])){
					temp.push_back(Mlv[posm]);
					posm++;
					break;
				}
				else if((pp!=mulpos[tt])&&(tt==(ed-3))){
					temp.push_back(Lev[posl]);
					posl++;
				}
			}
				
				}
	}
		rnd=rand() % (ed-1+el);
		posl=0;
				for(int pp=0;pp<rnd;pp++){
				Emv.push_back(temp[posl]);
				posl++;
				}
				Emv.push_back(b+k+h+1);
				
				for(int pp=rnd+1;pp<(ed-1+el);pp++){
				Emv.push_back(temp[posl]);
				posl++;
				}
			
			L.push_back(Emv);

}
		
// end of L segments;		
		
msize.resize(0);

for(int i=0;i<hz.size();i++){
	if(hz[i]!=0){
		for(int j=0;j<hz[i];j++){
			msize.push_back(i+1);
		}
	}
}

vector<int> Hsiz;

Hsiz=generateRandom(msize);
posl=k;


for(int i=0;i<Hsiz.size();i++){
	temp.resize(0);
	if((i+1)==posroot){
			temp.push_back(b+k+1); 
		}
	for(int j=0;j<Hsiz[i];j++){
		 
		temp.push_back(YY[posl]); //reuse the remaining part of YY as a random permutation of Z;
		posl++;
	}
	S.push_back(temp);
}

if((ed>2)||((ed==2)&&(el>0))){

temp.resize(0);
temp={b+k+1};
S.push_back(temp);
}
//end of S segments;



vector<int>  rootindex, tempL, tempS;
list<vector <int> > Sdot;


rootindex=pairroot( S,  L,  b,  k, h, ed, el);

for(int i=0; i<rootindex.size(); i++)
cout<< rootindex[i] <<'\n';

if((ed>2)||((ed==2)&&(el>0))){

for(int i=1; i<h+2; i++){
	temp.resize(0);
	tempL.resize(0);
	tempS.resize(0);
	tempS=git(S,i);
	tempL=git(L,i);
	
	for(int kk=0; kk<tempL.size(); kk++){
		cout << tempL[kk]<<'\n';
	}
	
	temp=rake2dot(tempS, tempL, b,k,h);
	
	for(int i=0; i< temp.size(); i++)
	cout << temp[i] << '\n';
	
	Sdot.push_back(temp);
}
}

else{

for(int i=1; i<h+1; i++){
	temp.resize(0);
	tempL.resize(0);
	tempS.resize(0);
	tempS=git(S,i);
	tempL=git(L,i);
	
	for(int kk=0; kk<tempL.size(); kk++){
		cout << tempL[kk]<<'\n';
	}
	
	temp=rake2dot(tempS, tempL, b,k,h);
	
	for(int i=0; i< temp.size(); i++)
	cout << temp[i] << '\n';
	
	Sdot.push_back(temp);
}
}

vector<int> dotbracket;


dotbracket=con_nest(Sdot, rootindex, b, k, h, posroot);

for(int i=0; i<dotbracket.size(); i++)
cout<< dotbracket[i];


cout<<'\n';
cout<<dotbracket.size()<<'\n';
cout<< b <<'\n';
cout<< k << '\n';


temp.resize(0);

for(int i=0; i<dotbracket.size(); i++)
temp.push_back(dotbracket[i]);

outdot(temp, netname);

/*
	ofstream stru;
	stru.open(netname+"_structure.txt");
	
	for(int i=0; i< temp.size(); i++){
	if(temp[i]==-1){
		stru<<'{';
	}
	
	else if(temp[i]==1){
		stru<<'}';
	}
	
	else if(temp[i]==0){
		stru<<'.';
	}
	
}

stru.close();

*/
}

return 0;

}
//end of main block  







