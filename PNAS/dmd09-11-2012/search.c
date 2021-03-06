#include <math.h>
#include <stdio.h>
#include <float.h>
#include "bcp.h"
#include "stdlib.h"
#include "search.h"
#include "controls.h"
#include "rng.h"
#define ERR1 (-1) 
#define ERR2 (-1) 
static tsearch search;
static atom * a;
static well_type ** ecoll;
static FILE * fp;

alist * treeadd(alist * p, alist * npt)
{  
  //search.b_d++;
  if(p==NULL)
      p=npt;
  else{
    npt->parent=p;
    if (npt->t<p->t)
      p->left=treeadd(p->left,npt);
    else if (npt->t>p->t)
      p->right=treeadd(p->right,npt);
    else if(rng()<0.5)
      p->left=treeadd(p->left,npt);
    else
      p->right=treeadd(p->right,npt);
  }
  return p;
}

alist * treeleft0(alist *p)
{
  //search.l_d++;
  if(p->left!=NULL)p=treeleft0(p->left);
  return p;
}


alist * treeleft(alist *p)
{
  if(p->left!=NULL)p=treeleft(p->left);
  return p;
}

alist * treeright(alist *p)
{
  if(p->right!=NULL)p=treeright(p->right);
  return p;
}

alist * treekill(alist * root,alist * npt)
{
  int branch=0;
  alist * parent=npt->parent;
  alist * left=npt->left;
  alist * right=npt->right;
  alist *root1=root;
  if(parent!=NULL)
    {
      if(parent->left==npt)
	branch=-1;
      else 
	branch=1;
    }

  if((left==NULL)&&(right==NULL))
    {
      if(branch>0)parent->right=NULL;
      else if(branch<0)parent->left=NULL;
      else root1=NULL;
    }
  else if (left==NULL) 
    {
      if(branch>0)parent->right=right;
      else if(branch<0)parent->left=right;
      else root1=right;
      right->parent=parent;
    }
  else if (right==NULL)
    {
      if(branch>0)parent->right=left;
      else if(branch<0)parent->left=left;
      else root1=left;
      left->parent=parent;
    }
  else if(rng()<0.5)
    {
      alist * right1=treeright(left);
      if(right1==left)
	{
	  if(branch>0)parent->right=left;
	  else if(branch<0)parent->left=left;
	  else root1=left;
	  left->parent=parent;
	  left->right=right;
	  right->parent=left;
	}
      else
	{
	  alist * parent1=right1->parent;
          alist * left1=right1->left;
	  parent1->right=left1;
	  if(left1!=NULL)left1->parent=parent1; 
	  if(branch>0)parent->right=right1;
	  else if(branch<0)parent->left=right1;
	  else root1=right1;
	  right1->parent=parent;
	  right1->right=right;
          right1->left=left;
          left->parent=right1;
          right->parent=right1; 
	}
    }
  else
   {
      alist * left1=treeleft(right);
      if(left1==right)
	{
	  if(branch>0)parent->right=right;
	  else if(branch<0)parent->left=right;
	  else root1=right;
	  right->parent=parent;
	  right->left=left;
          left->parent=right;
	}
      else
	{
	  alist * parent1=left1->parent;
          alist * right1=left1->right;
	  parent1->left=right1;
	  if(right1!=NULL)right1->parent=parent1; 
	  if(branch>0)parent->right=left1;
	  else if(branch<0)parent->left=left1;
	  else root1=left1;
	  left1->parent=parent;
	  left1->right=right;
          left1->left=left;
          left->parent=left1;
          right->parent=left1; 
	}
    }
  return root1;
}

void initsearch(void)
{ 
  alist * inpt;
  int i;
  int n=search.n;
  alist * tc=search.storage;
  alist ** ptc=search.begin;
  nlist * nc=search.neib_storage;
  nlist ** ntc=search.neib_begin;
  double eps,maxrb; 


  /*  fp=fopen("search_db","w");*/
  for (i=0;i<n;i++)
    {
    search.collp[i]=-1;
    search.collq[i]=-1;
    search.cell_storage[i].n=i;
    search.cell_storage[i].next=NULL;
    }
  search.np=0;
  search.nq=0;

  for(i=0;i<search.maxfree;i++)ptc[i]=&tc[i];
  ptc[search.maxfree]=NULL;
  search.free=search.begin;
for(i=0;i<search.maxfree;i++)ntc[i]=&nc[i];
  ntc[search.maxfree]=NULL;
  search.neib_free=search.neib_begin;
  //search.t_l=0;
  //search.b_d=0;
  //search.b_n=0;
  //search.l_n=0;
  //search.l_d=0;
  for(i=0;i<n;i++)
    {
      search.atom_collisionp[i]=NULL;
      search.atom_collisionq[i]=NULL;
    }
  for(i=0;i<=search.z;i++)
    search.cell_atoms[i]=NULL;
  for(i=0;i<n;i++)
    search.neib[i]=NULL;
  search.root=NULL;
}


int realloc_search (void)
{
  int i,j,k,level;
  int maxadd;
  dimensions * bound=get_bounds();
  if ((search.x==bound[0].period)&&(search.y==search.x*bound[1].period)
      &&(search.z==search.y*bound[2].period))return 1;
  search.x=bound[0].period;
  search.y=search.x*bound[1].period;
  if(search.z==search.y*bound[2].period)return 1;
  search.z=search.y*bound[2].period;
  maxadd=search.z;


if(maxadd>search.maxalloc)
  {
    free(search.cell_atoms);
    if(!(search.cell_atoms=(clist **)malloc((maxadd+1)*sizeof(clist *))))return 0;
  }
  return 1;
}

int allocsearch(int n)
{
 int i,j,k,level;
 int maxadd;
 dimensions * bound=get_bounds();
 a=get_atom();
 ecoll=get_ecoll();
 search.x=bound[0].period;
 search.y=search.x*bound[1].period;
 search.z=search.y*bound[2].period;
 maxadd=search.z;
 search.maxalloc=search.z;
 search.n=n;
 search.nhalf =n>>1;     
 if(maxadd>MAXADD)return 0;
 if(!(search.collp=(int *)malloc(n*sizeof(int))))return 0;
 if(!(search.collq=(int *)malloc(n*sizeof(int))))return 0;
 if(!(search.atomp=(int *)malloc(n*sizeof(int))))return 0;    
 if(!(search.atomq=(int *)malloc(n*sizeof(int))))return 0;
 if(!(search.fixed_atoms=(fatom *)malloc(n*sizeof(fatom))))return 0;

 
 if(!(search.cell_atoms=(clist **)malloc((maxadd+1)*sizeof(clist *))))return 0;

 if(!(search.cell_storage=(clist *)malloc(n*sizeof(clist))))return 0;

 search.maxfree=NFREE;
 if(!(search.storage =(alist *)malloc(search.maxfree*sizeof(alist))))return 0;
 if(!(search.begin=(alist **)malloc((search.maxfree+1)*sizeof(alist *))))return 0;
 if(!(search.neib=(nlist **)malloc((search.n)*sizeof(clist *))))return 0;
 if(!(search.neib_storage =(nlist *)malloc(search.maxfree*sizeof(nlist))))return 0;
 if(!(search.neib_begin=(nlist **)malloc((search.maxfree+1)*sizeof(nlist *))))return 0;

 if(!(search.atom_collisionp=(alist **)malloc((n)*sizeof(alist *))))return 0;
 if(!(search.atom_collisionq=(alist **)malloc((n)*sizeof(alist *))))return 0;
 
return 1;
}


/* used in change_neib*/
int find_collisions(int p1)
{
  int np=0;
  alist * pt=search.left;
  alist *ptp= pt->pprev;
  alist *ptn=pt->pnext;
  search.root=treekill(search.root,pt);
  *(--search.free)=pt;
  //search.t_l--;
  if(ptn!=NULL)ptn->pprev=ptp;  
  if(ptp!=NULL)ptp->pnext=ptn; 
  else
    search.atom_collisionp[p1]=ptn;
  pt=search.atom_collisionp[p1];
  while(pt)
    {
      int q=pt->q;
      if(search.collp[q]==-1)search.atomp[np++]=q;
      search.collp[q]=-2;
      pt=pt->pnext;
    }
  pt=search.atom_collisionq[p1];
  while(pt)
    {
      int q=pt->p;
      if(search.collp[q]==-1)search.atomp[np++]=q;
      search.collp[q]=-2;
      pt=pt->qnext;
    }
  return np;
}

/* used in squeeze_table */ 
int kill_collisions(int p1)
{
  alist ** free=search.free;
  alist * pt=search.atom_collisionp[p1];
  while(pt)
    {
      alist * qprev=pt->qprev;
      int q=pt->q;
      if(q<search.n)
	{
	  if(qprev==NULL)
	    {
	      search.atom_collisionq[q]=pt->qnext;
	      if(pt->qnext)
		search.atom_collisionq[q]->qprev=NULL;
	    }
	  else
	    { 
	      qprev->qnext=pt->qnext;
	      if(pt->qnext)
		pt->qnext->qprev=qprev;
	    }
	}
      search.root=treekill(search.root,pt);
      *(--free)=pt;
      //search.t_l--;
      pt=pt->pnext;
      if(pt==NULL)break;
      pt->pprev=NULL;
    }
  search.atom_collisionp[p1]=pt;
  pt=search.atom_collisionq[p1];
  while(pt)
    {
      alist * pprev=pt->pprev;
      int p=pt->p;
 	  if(pprev==NULL)
	    {
	      search.atom_collisionp[p]=pt->pnext;
	      search.atom_collisionp[p]->pprev=NULL;
	    }
	  else
	    { 
	      pprev->pnext=pt->pnext;
	      if(pt->pnext)
		pt->pnext->pprev=pprev;
	    }
     search.root=treekill(search.root,pt);
      *(--free)=pt;
      //search.t_l--;
      pt=pt->qnext;
      if(pt==NULL)break;
      pt->qprev=NULL;
    }
  search.atom_collisionq[p1]=pt;   
  search.free=free;
}


double get_coll_time(int p1,int q1)
{
  alist * pt=search.atom_collisionp[p1];
  while(pt)
    {
      if(pt->q==q1)return pt->t; 
      pt=pt->pnext;
    }
 pt=search.atom_collisionq[p1];
  while(pt)
    {
      if(pt->p==q1)return pt->t; 
      pt=pt->qnext;
    }
  return -100.0;
}




/* used in new_loc in bcp.c */
void change_cell(int p1, size_al address1)
{ 
  alist * pt=search.left;
  alist *ptp= pt->pprev;
  alist *ptn=pt->pnext;
  clist * inpt;
  int address=a[p1].add;
  clist * old_rec=search.cell_atoms[address];
  if(old_rec->n==p1)
    search.cell_atoms[address]=search.cell_atoms[address]->next;
  else
    {
      clist *  prev_rec=old_rec;
      old_rec=prev_rec->next;
      while(old_rec->n!=p1)
	{
	  prev_rec=old_rec;

	  old_rec=prev_rec->next;
	}
      prev_rec->next=old_rec->next;
    }
  old_rec->next=search.cell_atoms[address1];
  search.cell_atoms[address1]=old_rec;
  a[p1].add=address1;

  search.root=treekill(search.root,pt);
  *(--search.free)=pt;
  //search.t_l--;
  if(ptn!=NULL)ptn->pprev=ptp;  
  if(ptp!=NULL)ptp->pnext=ptn; 
  else
    search.atom_collisionp[p1]=ptn;
  search.np=0;
  search.nq=0;
}

void change_neib(int p1)
{ 
  clist * inpt;
  int nq;
  int p2;
  int i,j,k,i1,j1,k1,i2,j2,k2,address=a[p1].add;
  int addressz,addressy;
  int ct;
  alist * pt=search.left;
  alist *ptp= pt->pprev;
  alist *ptn=pt->pnext;
  search.root=treekill(search.root,pt);
  *(--search.free)=pt;
  if(ptn!=NULL)ptn->pprev=ptp;  
  if(ptp!=NULL)ptp->pnext=ptn; 
  else
    search.atom_collisionp[p1]=ptn;
  search.np=0;
  search.nq=0;
  if(p1==ERR1)
    {
      int c1=get_coll_type(ERR1,ERR2);
      int c2=get_coll_type(ERR2,ERR1);
      printf("in %d %d\n",c1,c2);
    }

  nq=kill_neib_list(p1);

  i1=a[p1].i.x.i-1;
  j1=(a[p1].i.y.i-1)*search.x;
  i2=i1+2;
  j2=j1+(search.x<<1);
  if(search.z==search.y)
  {k1=a[p1].i.z.i*search.y;k2=k1;}
  else
  {k1=(a[p1].i.z.i-1)*search.y;k2=k1+(search.y<<1);}
  
  for(k=k1;k<=k2;k+=search.y)
    { 
      addressz=k;
      if(addressz<0)addressz+=search.z;
      if(addressz==search.z)addressz=0;
      for(j=j1;j<=j2;j+=search.x)
	{ 
	  addressy=j; 
	  if(addressy<0)addressy+=search.y;
	  if(addressy==search.y)addressy=0;
	  addressy+=addressz;
	  for(i=i1;i<=i2;i++)
	    { 
	      address=i; 
	      if(address<0)address+=search.x;
	      if(address==search.x)address=0;
	      address+=addressy;
	      for(inpt=search.cell_atoms[address];inpt;inpt=inpt->next)
		{
		  int p2=inpt->n;
		  if(p1!=p2)
		    if(atom_dist2(p1,p2)<search.rneib)
		      {
			if(search.collq[p2]<-1)
			    search.collq[p2]=-2-search.collq[p2];			  
			else
			  { 
			    search.atomq[nq++]=p2;
			    search.collq[p2]=ecoll[a[p1].c][a[p2].c];
			  }
		      }
		}
	    } 
	}
    }
  for(i=0;i<nq;i++)
    {
      p2=search.atomq[i];
      ct=search.collq[p2];
      if(ct>-1)
	{ 
	  add_neib(p1,p2,ct);
	  if(get_coll_type(p2,p1)<=-1)
	    add_neib(p2,p1,ct);
	}
      else if(ct<-1)
	kill_neib(p2,p1);
      search.collq[p2]=-1;
    }
  if(p1==ERR1)
    {
      int c1=get_coll_type(ERR1,ERR2);
      int c2=get_coll_type(ERR2,ERR1);
      printf("out %d %d\n",c1,c2);
    }

}

int kill_neib_list(int i)
{
  int nq=0;
  nlist * pt=search.neib[i];
  while (pt)
    {
      int q=pt->n;
      if(search.collq[q]==-1)
	{
	  search.atomq[nq++]=q;  
	  search.collq[q]=-2-(pt->c);
	}  
      *(--search.neib_free)=pt;
      pt=pt->next;
    }
  search.neib[i]=NULL;
  return nq;
}

int get_coll_type(int i,int j)
{
  nlist * pt=search.neib[i];
  while (pt)
    {
      if(pt->n==j)
	return pt->c;
      pt=pt->next;
    }
  return -1;
}

int set_coll_type(int i,int j,int ct)
{
  nlist * pt=search.neib[i];
  while (pt)
    {
      if(pt->n==j)
	{
	  pt->c=ct;
	  return pt->c;
	}
      pt=pt->next;
    }
  return -1;
}

int kill_neib(int i, int j)
{
  nlist * pt=search.neib[i];
  if(pt)
    {
      nlist * pt1=pt;
      if(pt->n==j)
	{
	  search.neib[i]=pt->next;
	  *(--search.neib_free)=pt;
	  return 1;
	}
      pt=pt->next;
      while (pt)
	{
	  if(pt->n==j)
	    {
	      pt1->next=pt->next;
	      *(--search.neib_free)=pt;
	      return 1;
	    }
	  pt1=pt;
	  pt=pt->next;
	}
    }
  return 0;
}


/* find atoms in ajacent cells to a given atom; 
   returns number of such atoms */
int find_neighbors(int p, int * neib)
{ 
 clist * inpt;
 int nn=0;
 int i,j,k,i1,j1,k1,i2,j2,k2,address;
 int addressz,addressy;
  
  i1=a[p].i.x.i-1;
  j1=(a[p].i.y.i-1)*search.x;
  i2=i1+2;
  j2=j1+(search.x<<1);
  if(search.z==search.y)
  {k1=a[p].i.z.i*search.y;k2=k1;}
  else
  {k1=(a[p].i.z.i-1)*search.y;k2=k1+(search.y<<1);}
  
  for(k=k1;k<=k2;k+=search.y)
    { 
      addressz=k;
      if(addressz<0)addressz+=search.z;
      if(addressz==search.z)addressz=0;
      for(j=j1;j<=j2;j+=search.x)
	{ 
	  addressy=j; 
	  if(addressy<0)addressy+=search.y;
	  if(addressy==search.y)addressy=0;
	  addressy+=addressz;
	  for(i=i1;i<=i2;i++)
	    { 
	      address=i; 
	      if(address<0)address+=search.x;
	      if(address==search.x)address=0;
	      address+=addressy;
	      for(inpt=search.cell_atoms[address];inpt;inpt=inpt->next)
		    neib[nn++]=inpt->n;
	    } 
	}  
    }
 return nn;
}  
/* list atoms in a certain cell, returns number of atoms in a cell */
int list_atoms(size_al address,int * atomx)
{ 
	int i=0;
	clist *inpt1;
	for( inpt1=search.cell_atoms[address];inpt1;inpt1=inpt1->next)
 	{ atomx[i++]=inpt1->n;}
 return i;	
}

int bond(int i, int j)
{   
    nlist *inpt1;
    for( inpt1=search.neib[i];inpt1;inpt1=inpt1->next)
      if(j==inpt1->n) return inpt1->c;
    return -ecoll[a[i].c][a[j].c];
}

void add_atom2cell(int i)
{
  size_al address=a[i].add;
  clist * pt=search.cell_atoms[address];
  search.cell_atoms[address]=&search.cell_storage[i];
  search.cell_atoms[address]->next=pt;
}


int init_tables(void)
{
  int i0,j0,k0,address,ix,iy,iz,level;
  int i,j,k,i1,j1,k1,i2,j2,k2;
  int maxadd=search.z;
  size_al address1;
  size_al address2;
  int n=search.n-1;
  int n1=n+1;
  double t;
  int p,q,ct;
  int addressz,addressy;
  fatom * fa=search.fixed_atoms;
  nlist * pt;
  initsearch();
  for(i=0;i<=n;i++)
    {
      add_atom2cell(i);
      q=twall(i,&t);
      bubble(t,i,q);
      /* initializing original positions of all the atoms for neiborlists */ 
      fa[i].x=a[i].r.x;
      fa[i].y=a[i].r.y;
      fa[i].z=a[i].r.z;
      t=tneib(i);
      bubble(t,i,search.n+3);
    }
  
  
  a[n1].c=0;
  for (k0=0;k0<maxadd;k0++)
    if(search.np=list_atoms(k0,search.atomp))
      {
	i0=search.atomp[0];
	i1=a[i0].i.x.i-1;
	i2=i1+2;
	j1=(a[i0].i.y.i-1)*search.x;
	j2=j1+(search.x<<1);
	if(search.z==search.y)
	  {k1=a[i0].i.z.i*search.y;k2=k1;}
	else
	  {k1=(a[i0].i.z.i-1)*search.y;k2=k1+(search.y<<1);}
	for(k=k1;k<=k2;k+=search.y)
	  { 
	    addressz=k;
	    if(addressz<0)addressz+=search.z;
	    if(addressz==search.z)addressz=0;
	    for(j=j1;j<=j2;j+=search.x)
	      { 
		addressy=j; 
		if(addressy<0)addressy+=search.y;
		if(addressy==search.y)addressy=0;
		addressy+=addressz;
		for(i=i1;i<=i2;i++)
		  { 
		    int iq;
		    int ip;
		    address=i; 
		    if(address<0)address+=search.x;
		    if(address==search.x)address=0;
		    address+=addressy;
		    search.nq=list_atoms(address,search.atomq);
		    for (iq=0;iq<search.nq;iq++)
		      for (ip=0;ip<search.np;ip++)
			{
			  i0=search.atomp[ip]; 
			  j0=search.atomq[iq];
			  if(i0<j0)			
			    if(atom_dist2(i0,j0)<search.rneib)
			      {
				int ct=collision_type(i0,j0);
				add_neib(i0,j0,ct);
				add_neib(j0,i0,ct);    
			    }			
		       }
		  }
	      }
	  }
      }
 
  for(i0=0;i0<=n;i0++)
    {
      pt=search.neib[i0];
      while(pt)
	{
	  j0=pt->n;	      
	  if(i0<j0)
	    {
	      ct=pt->c;
	      if(ct<0) return (int) ct;
	      add_potential(ct);
	      if(tball(i0,j0,ct,&t))
		bubble(t,i0,j0);
	    }
	  pt=pt->next;
	}
    }
  return 1;  
}

int bubble (double t, int p1, int q1)
     //bubble is no longer doing bubble sort, it does treesort
{
  int p,q;
  alist *npt,*pt,*inpt;
  if(t<0)
    printf("bubble %le %d %d\n",t,p1,q1);
  if(q1<p1)
    {
      p=q1;
      q=p1;
    }
  else
    {
      p=p1;
      q=q1;
    }
  if((q<search.n)&&(q-p>search.nhalf))
    {
      int r=q;
      q=p;
      p=r;
    }

  if (p>=0)
    {
      npt=*(search.free);
      if(!npt){writetext(get_text_name());exit(0);}
      search.free++;
      //search.t_l++;
      //search.b_n++;  
      npt->t=t;
      npt->p=p;
      npt->q=q;
      npt->parent=NULL;
      npt->left=NULL;
      npt->right=NULL;

      search.root=treeadd(search.root,npt);

      if(q<search.n)
	{
	  inpt=search.atom_collisionq[q];	
	  npt->qnext=inpt;
          if(inpt)  
	    inpt->qprev=npt;          
	  npt->qprev=NULL;
	  search.atom_collisionq[q]=npt;
	}
      inpt=search.atom_collisionp[p];
      npt->pnext=inpt;
      if(inpt)
	inpt->pprev=npt;          
      npt->pprev=NULL;
      search.atom_collisionp[p]=npt;
      return p;
    }
  return -1;
}


int add_neib (int p, int q, int c)
{
  nlist *npt,*pt,*inpt;
  npt=*(search.neib_free);
  if(!npt){writetext(get_text_name());exit(0);}
  search.neib_free++;
  npt->n=q;
  npt->c=c;
  npt->next=search.neib[p];
  search.neib[p]=npt;
  return 1;
}



int get_free(void)
{
  return (int)(search.free-search.begin);
}
int get_maxfree(void)
{
  return search.maxfree;
}
void set_maxfree(int a)
{
  search.maxfree=a;
}

fatom * get_fixed_atom(int i){
  return (search.fixed_atoms)+i;
}


void update_table(int p1, int q1, int ct1)
{
  int i,p2,ct;
  int add;
  double t;
  int q;
  nlist * npt;
  if(q1<search.n+3)
    {
      q=twall(p1,&t);
      bubble(t,p1,q);
    }
  if(q1==search.n+3)
    {
      search.fixed_atoms[p1].x=a[p1].r.x;
      search.fixed_atoms[p1].y=a[p1].r.y;
      search.fixed_atoms[p1].z=a[p1].r.z;
      t=tneib(p1);
      bubble(t,p1,q1);
     npt=search.neib[p1];
      while(npt){
	p2=npt->n;
	ct=npt->c;
	if(tball(p2,p1,ct,&t)) 
	  bubble(t,p1,p2);
	npt=npt->next;
      } 
   }
  if (q1<search.n)
    {
      set_coll_type(p1,q1,ct1);
      set_coll_type(q1,p1,ct1);
      t=tneib(p1);
      bubble(t,p1,search.n+3);
      q=twall(q1,&t);
      bubble(t,q1,q);
      t=tneib(q1);
      bubble(t,q1,search.n+3);
      npt=search.neib[p1];
      while(npt){
	p2=npt->n;
	ct=npt->c;
	if(tball(p2,p1,ct,&t)) 
	  bubble(t,p1,p2);
	npt=npt->next;
      }
      npt=search.neib[q1];
      while(npt){
	p2=npt->n;
	if(p2!=p1)
	  {
	    ct=npt->c;
	    if(tball(p2,q1,ct,&t)) 
	      bubble(t,q1,p2);
	  }
	npt=npt->next;
      }
    }
}




int squeeze_table(int * p1, int * q1, double * timea)
{ 
  alist *inpt=treeleft0(search.root);
  search.left=inpt;
  //  search.l_n++;
  //printf("%20.1l6f %d %d %d\n",inpt->t,inpt->p,inpt->q,inpt->c);
  *timea=inpt->t;
  *p1=inpt->p;
  *q1=inpt->q;
  if(*timea<0)
  //   printf("%d %d %le\n",*p1,*q1,*timea);
   if (((*p1==ERR1)&&(*q1==ERR2))||((*q1==ERR1)&&(*p1==ERR2))) 
   {
     int c1=get_coll_type(ERR1,ERR2);
     int c2=get_coll_type(ERR2,ERR1);
     double t1=get_coll_time(ERR1,ERR2);
     double t2=get_coll_time(ERR2,ERR1);
     double w1=a[ERR1].w;
     double w2=a[ERR2].w;
     double wt1=get_coll_time(ERR1,1731);
     double wt2=get_coll_time(ERR2,1731);
     printf("%d %d %le %d %d %le %le %le %le %le %le\n",*p1,*q1,*timea,c1,c2,t1,t2,w1,w2,wt1,wt2);
     } 
if((*q1)>=search.n)
    return -1;
  /*kill_collisions works only if it is an atomic collision */
   else
    {
      kill_collisions(*p1);
      kill_collisions(*q1);
    }
  return get_coll_type(*p1,*q1);
}


int * get_collp(void){return search.collp;}
int * get_collq(void){return search.collq;}
int * get_atomp(void){return search.atomp;}
int * get_atomq(void){return search.atomq;}



int get_nq(int i)
{
  int nq=0;
  nlist * pt=search.neib[i];
  while (pt)
    {
      int q=pt->n;
      if(search.collq[q]==-1)
	{
	  search.atomq[nq++]=q;  
	  search.collq[q]=pt->c;
	}  
    }
  search.nq=nq;
  return nq;
}

int get_np(int i)
{
  int nq=0;
  nlist * pt=search.neib[i];
  while (pt)
    {
      int q=pt->n;
      if(search.collp[q]==-1)
	{
	  search.atomp[nq++]=q;  
	  search.collp[q]=pt->c;
	}  
    }
  search.np=nq;
  return nq;
}
void clean_neib(void)
{
  int i;
  int np=search.np;
  for(i=0;i<np;i++)
    search.collp[search.atomp[i]]=-1;
  np=search.nq;
  for(i=0;i<np;i++)
    search.collq[search.atomq[i]]=-1;
  search.np=0;
  search.nq=0;
}

void update_neib(int i1, int i2)
{
  int np,j,ct;
  nlist * npt=search.neib[i1];
  int i=0;
  while(npt)
    {
      j=search.atomp[i];
      ct=search.collp[j];
      npt->c=ct;
      set_coll_type(j,i1,ct);
      i++;
      npt=npt->next;
    }
  if(i!=search.np)
    printf("error in update neibors p %d %d\n",i,search.np);

  npt=search.neib[i2];
  i=0;
  while(npt)
    {
      j=search.atomq[i];
      ct=search.collq[j];
      npt->c=ct;
      set_coll_type(j,i2,ct);
      i++;
      npt=npt->next;
    }
  if(i!=search.nq)
    printf("error in update neibors q %d %d\n",i,search.nq);
  clean_neib();
}




int pairs(int (*do_something)(int,int,int))
{ 
  int n_bonds=0;
  int i,j;
  nlist *inpt1;
  for(i=0;i<search.n;i++)
    for(inpt1=search.neib[i];inpt1;inpt1=inpt1->next)
      { 
	j=inpt1->n;
	n_bonds+=do_something(i,j,inpt1->c);
      }
  return n_bonds;
}  

void print_collisions(int p1,int q1)
{
  int np=0;
  alist * pt=search.atom_collisionp[p1];
      fprintf(fp,"collision %lf %d %d %d\n",get_time(),get_ll(),p1,q1);
      while(pt)
	{
	  fprintf(fp,"%10.6lf %3d %3d %2d\n",pt->t,pt->p,pt->q,get_coll_type(pt->p,pt->q));
	  pt=pt->pnext;
	}
      pt=search.atom_collisionq[p1];
      fprintf(fp,"--------\n");
      while(pt)
	{
	  fprintf(fp,"%10.6lf %3d %3d %2d\n",pt->t,pt->p,pt->q,get_coll_type(pt->p,pt->q));
	  pt=pt->qnext;
	}
      fflush(fp);
}

int find_shell(int p1, int ct, int * neib)
{ 
  int nn=0;
  nlist * pt=search.neib[p1];
  while(pt)
    {
      if(pt->c==ct)
	{
	  neib[nn]=pt->n;
	  nn++;
	}
      pt=pt->next;
    }
  return nn;
}

void set_rneib(double eps1, double maxrb1)
{
  double eps=eps1;
  double maxrb=maxrb1;
  if(eps<=0)eps=maxrb/6;
  search.eps2=eps*eps;
  maxrb+=3*eps;
  search.rneib=maxrb*maxrb;  
  printf("maxrb=%lf eps=%lf\n", maxrb, eps);
}
double get_eps2(void){return search.eps2;}

