#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include"ran2.c"
#define L 2000
long sd=-937176787;
main()
{
  int i,j,count2,iw,count,n,NRUN=3,t,TMAX=10000,ii,jj,ni,nj,id,xl,yl,M=50,tt,ratio,RELAX=50000,counti,countj,avc_dist,up,down;
  int super_spin[L],zero,s_up,s_down,s_zero,avc_a,diff_spin,avc_opdist=0;
  float stat[L],sum,sum_sq,sum_forth,sum_tot,BQ,BQ_tot,sum1,sum2,sum3,sum4,suspt,r,p,mu,op[TMAX],tot_sum,s_op[TMAX],time_avg[TMAX],tot_flip,tot_zero,mag,dist_down[L],dist_up[L],dist_s_up[L],dist_zero[L],op_pvt[L], dist_s_down[L], dist_s_zero[L],tot_s_zero,flag1,dist_op[L],flag,store,delta,k;
  //FILE *fpt1,*fpt2;

FILE *fp;

  delta=0.60;
  k=1.5;

 //   p=1.0;

  for(i=0;i<L;i++) 
  {
    stat[i]=0; 
    dist_up[i]=0; dist_down[i]=0; dist_zero[i]=0;dist_op[i]=0;
   /* if(i<M)*/ {super_spin[i]=0; dist_s_up[i]=0; dist_s_down[i]=0; dist_s_zero[i]=0;}
 

  }

//  for(M=1;M<=L/2;M=M*2)
//for(delta=0;delta<1.;delta=delta+0.01)
{
 for(p=0;p<0.4;p=p+0.01)
 {
  M=50;
  
  avc_dist=0;
  tot_s_zero=0;
  ratio=L/(M);
  //for(i=0;i<M;i++){for(j=0;j<M;j++) {super_spin[i][j]=0;}}

  for(i=0;i<TMAX;i++) {op[i]=0; s_op[i]=0; time_avg[i]=0;}

/////////////relaxation//////////////////////////
     for(i=0;i<L;i++)
     {
       
          if(ran2(&sd)<1.) {stat[i]=1;}
          else {stat[i]=-1;}
        op_pvt[i]=2*ran2(&sd)-1;
     } 

    for(t=0;t<RELAX;t++)
    {
        for(j=0;j<L;j++)
         {
            ii=ran2(&sd)*L; jj=ran2(&sd)*L;
 

            if(ran2(&sd)<p) {mu=-1;} else {mu=+1;}
            store=stat[ii];
            stat[ii]=stat[ii]+mu*stat[jj];
            if(stat[ii]<-1) {stat[ii]=-1;}
            if(stat[ii]>1) {stat[ii]=1;}
            op_pvt[ii]=op_pvt[ii]+k*(stat[ii]-store);
            if(op_pvt[ii]<-1) {op_pvt[ii]=-1;}
            if(op_pvt[ii]>1) {op_pvt[ii]=1;}

            if(fabs(op_pvt[ii]-stat[ii])>delta) 
            {
               if(op_pvt[ii]<0) {stat[ii]=-1;}
               if(op_pvt[ii]>0) {stat[ii]=1;}
            } 
            

         }
    }

/////////////////////////////////////////////////

  zero=0; tot_flip=0; flag1=0; tot_zero=0; BQ_tot=0.0; sum_tot=0.0;

  for(n=0;n<NRUN;n++)
  {

     tt=0;
     for(t=0;t<TMAX;t++)
     {
         for(j=0;j<L;j++)
         {
              ii=ran2(&sd)*L; jj=ran2(&sd)*L;
 

            if(ran2(&sd)<p) {mu=-1;} else {mu=+1;}
            store=stat[ii];
            stat[ii]=stat[ii]+mu*stat[jj];
            if(stat[ii]<-1) {stat[ii]=-1;}
            if(stat[ii]>1) {stat[ii]=1;}
            op_pvt[ii]=op_pvt[ii]+k*(stat[ii]-store);
            if(op_pvt[ii]<-1) {op_pvt[ii]=-1;}
            if(op_pvt[ii]>1) {op_pvt[ii]=1;}

        //     printf("%f\n",op_pvt[ii]);

            if(fabs(op_pvt[ii]-stat[ii])>delta) 
            {
               if(op_pvt[ii]<0) {stat[ii]=-1;}
               if(op_pvt[ii]>0) {stat[ii]=1;}
            } 
         }

 
        ii=1; jj=0; 
         zero=0; counti=0; mag=0; flag=0; up=0; down=0; s_up=0; s_down=0; s_zero=0;
         for(i=0;i<L;i++)
         {
            mag+=stat[i];
            if(stat[i]==0) {zero++;}
            else if (stat[i]==1) {up++;}
            else {down++;}
  
         }
         counti=0; ii=0; avc_a=0;
         for(i=0;i<L;i++)
         {
            counti+=stat[i]; ii++; avc_a++;
            if(ii==ratio /*|| ii==(L-1)*/)
            {
              if(counti>0) {flag=flag+ran2(&sd); s_up++;}
              if(counti<0) {flag=flag-ran2(&sd); s_down++;}
              if(counti==0) {s_zero++;}
              counti=0; ii=0;
            //  if(avc_a!=ratio) {printf("%d\n",avc_a);}
              avc_a=0;
            //  printf("%f %d\n",stat[i],flag);
            }
         }
         tot_zero+=(float)zero/L;

         op[t]+=fabs(mag)/L;

   //     printf("%d %d %d\n",up,down,zero);

        dist_up[up]++;
        dist_down[down]++;
        dist_zero[zero]++;

        dist_s_up[s_up]++;
        dist_s_down[s_down]++;
        dist_s_zero[s_zero]++;       
        tot_s_zero+=(float)s_zero/M;
        avc_dist++;

         s_op[t]+=(float)flag/M;
      //  printf("%f %d\n",mag,flag);
         if(mag*flag<0) {flag1++; }
      //  if(mag>0 && flag>0) {flag1++;}
      //   flag1+=fabs(mag);
         diff_spin=fabs(mag);
         dist_op[diff_spin]++; avc_opdist++;
     }
      sum=0;
      sum_sq=0;
      sum_forth=0;
      for(ii=0;ii<TMAX;ii++) {sum+=op[ii];}
      for(ii=0;ii<TMAX;ii++) {sum_sq+=pow(op[ii],2);}
      for(ii=0;ii<TMAX;ii++) {sum_forth+=pow(op[ii],4);}
      
      sum=sum/TMAX;
      
      sum1=pow(sum,2);
      sum2=sum_sq/TMAX;
      sum3=pow(sum2,2);
      sum4=sum_forth/TMAX;
      BQ=1-((sum4)/(3*sum3));
      suspt=L*(sum2-sum1);

  }
//end of realisation

	fp=fopen("orderParameter_vs_p_k1.50_delta0.60_L2000.dat", "a");
	fprintf(fp,"%f\t%f\t%f\t%f\t%f\t%f\n",p,delta,sum/NRUN,tot_zero/(TMAX*NRUN),(float)flag1/(NRUN*TMAX),tot_s_zero/(NRUN*TMAX));
	fclose(fp);
	
	
 }
//end of p loop
///////////// printf("\n");
}
//end of delta loop
    
 // for(i=0;i<L;i++) {printf("%d %f %f %f %f\n",i,dist_zero[i]/avc_dist, dist_up[i]/avc_dist, dist_down[i]/avc_dist, dist_op[i]/avc_opdist);}

////  for(i=0;i<M;i++) {printf("%d %f %f %f %f\n",i,dist_s_zero[i]/avc_dist, dist_s_up[i]/avc_dist, dist_s_down[i]/avc_dist,tot_s_zero/avc_dist);}

  for(i=0;i<TMAX;i++)
  {
  //   fprintf(fpt1,"%d %f %f\n",i,op[i]/NRUN,s_op[i]/NRUN);
  //   fprintf(fpt2,"%d %f\n",i,time_avg[i]/NRUN);
  }


}
