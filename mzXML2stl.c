/*                                                                                                 */      /*  mzXML2stl - convert (region of) mzXML file to STL                                              */
/*                                                                                                 */
/* (c) Magnus Palmblad, The University of Reading, 2007-                                           */ 
/*                                                                                                 */
/*                                                                                                 */
/* This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;       */
/* without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.       */
/*                                                                                                 */
/* Contact information: magnus.palmblad@gmail.com                                                  */
/*                                                                                                 */
/*                                                                                                 */
/* compile with e.g. gcc -o mzXML2stl mzXML2stl.c base64.c ramp.c -I. -lgd -lm -lz                 */
/*                                                                                                 */
/*                                                                                                 */

#include <stdio.h>
#include <stdlib.h>  
#include <ctype.h>
#include <string.h>
#include <math.h>
#include "ramp.h"
#include "base64.h"



/* mzXML read functions adopted from Pep3D */

typedef struct {
  int size;
  double * xval;
  double * yval;
} spectStrct;

void freeSpectStrct(spectStrct spectrum)
{
  free(spectrum.xval);
  free(spectrum.yval);
  
  return;
}

ramp_fileoffset_t *pScanIndex;
int iLastScan;
struct ScanHeaderStruct scanHeader;
struct RunHeaderStruct runHeader;
ramp_fileoffset_t indexOffset;
spectStrct *spects;
double *wghs;
int spectNum;
int **z, **stl_height, **smooth, **smooth2, *new_spectrum;
long i;
RAMPREAL *pPeaks;
int n, n_MS_spectra=0;
RAMPFILE *pFI;

float SG[5][5] = {{-0.0743,0.0114,0.0400,0.0114,-0.0743},
		  {0.0114,0.0971,0.1257,0.0971,0.0114},
		  {0.0400,0.1257,0.1543,0.1257,0.0400},
		  {0.0114,0.0971,0.1257,0.0971,0.0114},
		  {-0.0743,0.0114,0.0400,0.0114,-0.0743}};

int getMsSpect(spectStrct *msSpect, RAMPFILE *pFI, int scanNum[2])
{
  void getCmbSpect(spectStrct *cmbSpect, int spectNum, spectStrct *spects, double *wghs);

  msSpect->size = -1;
  scanNum[0]=scanNum[0]>1?scanNum[0]:1;
  scanNum[1]=scanNum[1]<iLastScan?scanNum[1]:iLastScan;
  spectNum=scanNum[1]-scanNum[0]+1;
  if(spectNum < 1){
    printf("invalid scan number: %d-%d (full scan range: 1-%d)\n",scanNum[0], scanNum[1], iLastScan);
    fflush(stdout);
    free (pScanIndex);
    return -1;    
  }

  spects=(spectStrct *) calloc(spectNum, sizeof(spectStrct));
  spectNum=0;
  for(i=scanNum[0];i<=scanNum[1];++i) 
    {
      if((scanHeader.msLevel==1)&&(scanHeader.peaksCount>0)) /* MS ? */
	{                         
	  spects[spectNum].size=scanHeader.peaksCount;
	  spects[spectNum].xval=(double *) calloc(spects[spectNum].size, sizeof(double));
	  spects[spectNum].yval=(double *) calloc(spects[spectNum].size, sizeof(double));
	  
	  pPeaks=readPeaks(pFI,pScanIndex[i]);
	  
	  spects[spectNum].size=0;
	  n = 0;
	  while(pPeaks[n]!=-1)
	    {
	      spects[spectNum].xval[spects[spectNum].size]=pPeaks[n];
	      n++;
	      spects[spectNum].yval[spects[spectNum].size]=pPeaks[n];
	      n++;
	      ++(spects[spectNum].size);
	    }
	  free (pPeaks);
	  n_MS_spectra++;
	  if(spects[spectNum].size>0) ++spectNum; 
	  else freeSpectStrct(spects[spectNum]);
	}
      
    } 
  
  if(spectNum>0) 
    {
    wghs= (double *) calloc(spectNum, sizeof(double));
    for (i=0;i<spectNum;++i)
      wghs[i]=1.;
    getCmbSpect(msSpect, spectNum, spects, wghs);
    free(wghs);
  }
  else 
    {
      printf("cannot find an MS spectrum..."); fflush(stdout);
      for(i=0;i<spectNum;++i) freeSpectStrct(spects[i]);
      free(spects);
      return -1;
    }
  
  for(i=0;i<spectNum;++i) freeSpectStrct(spects[i]);
  free(spects);
  
  return 0;
}

void getCmbSpect(spectStrct *cmbSpect, int spectNum, spectStrct *spects, double *wghs)
{
  void copySpectStrct(spectStrct * tgtSpect, spectStrct srcSpect);
  spectStrct tmpSpect[2];
  int indx, indx1, indx2;
  double tmpWghs[2]={1.,1.};
  int i;
  
  if(spectNum<1) return;
  if(spectNum==1) 
    {
      copySpectStrct(cmbSpect,spects[0]);
      if(wghs[0]!=1.) for(i=0;i<cmbSpect->size;++i) cmbSpect->yval[i]*=wghs[0];
      return;
    } 
  
  return;
}

void copySpectStrct(spectStrct * tgtSpect, spectStrct srcSpect)
{
  int i;

  tgtSpect->size=srcSpect.size;
  tgtSpect->xval=(double *) calloc(tgtSpect->size, sizeof(double));
  tgtSpect->yval=(double *) calloc(tgtSpect->size, sizeof(double));

  for(i=0;i<tgtSpect->size;++i) 
    {
    tgtSpect->xval[i]=srcSpect.xval[i];
    tgtSpect->yval[i]=srcSpect.yval[i];
    }

  return;
}



/* main program starts here */

int main(int argc, char *argv[]) 
{
  FILE *inp, *outp;
  spectStrct mzXML_spectrum, temp_spectrum;
  RAMPFILE *mzXML_file;
  char *p, mzXML_filename[100], stl_filename[100], scan_range[100], mass_range[100], pixels[100], s[200];
  int range[2];
  ramp_fileoffset_t offset, *scan_index;
  struct ScanHeaderStruct scan_header;
  struct RunHeaderStruct run_header;

  int *scan_number;
  int XWIDTH,YWIDTH;
  
  int mzXML_MS_scan, scan, first_scan, last_scan, n_scans, x, y, c, lowest_mass, highest_mass, mz_range, horizontal_pixels, vertical_pixels, n_horizontal_labels, n_vertical_labels, horizontal_label_interval, vertical_label_interval, first_horizontal_label, first_vertical_label, MS_scans, noaxes;
  long i,j,k;
  float min_yval, max_yval, log_min_yval, log_max_yval, power, mz_per_pixel, scans_per_pixel, background, saturation, xscale, yscale, zscale, zbase, upper, lower, xmin, xmax, box_height, box_base;

 
  /* parsing command line parameters */
  
  if( (argc==2) && ( (strcmp(argv[1],"--help")==0) || (strcmp(argv[1],"-help")==0) || (strcmp(argv[1],"-h")==0)) ) /* want help? */
    {
      printf("mzXML2stl - (c) Magnus Palmblad 2007-\n\nusage: mzXML2stl -i<mzXML file> -o<STL filename> -s<first scan>,<last scan> -m<lowest mass>,<highest mass> [-c<color scheme> -p<horizontal pixels>,<vertical pixels> -B<background (absolute value)> -S<saturation (as fraction of maximum)>].\n\nfor more information, see http://www.ms-utils.org/mzXML2stl or e-mail magnus.palmblad@gmail.com\n");
      return 0;
    }
  
  if (argc<5 || argc>22) /* incorrect number of parameters */
    {
      printf("usage: mzXML2stl -i<mzXML file> -o<STL filename> -s<first scan>,<last scan> -m<lowest mass>,<highest mass> [-p<horizontal pixels>,<vertical pixels> -B<background (absolute value)> -S<saturation (as fraction of maximum)>] (type mzXML2stl --help for more information)\n");
      return -1;
    }
  
  
  /* set default values and parse flags */

  background=0.0; saturation=1.0; noaxes=1; /* linear, m/z on horizontal axis, noaxis yet */
  
  for(i=1;i<argc;i++) {
    if( (argv[i][0]=='-') && (argv[i][1]=='i') ) strcpy(mzXML_filename,&argv[strlen(argv[i])>2?i:i+1][strlen(argv[i])>2?2:0]);
    if( (argv[i][0]=='-') && (argv[i][1]=='o') ) strcpy(stl_filename,&argv[strlen(argv[i])>2?i:i+1][strlen(argv[i])>2?2:0]);
    if( (argv[i][0]=='-') && (argv[i][1]=='s') ) 
      {
	strcpy(scan_range,&argv[strlen(argv[i])>2?i:i+1][strlen(argv[i])>2?2:0]); p=strtok(scan_range,","); 
	first_scan=atoi(p); p=strtok('\0',","); last_scan=atoi(p);
      }
    if( (argv[i][0]=='-') && (argv[i][1]=='m') ) 
      {
	strcpy(mass_range,&argv[strlen(argv[i])>2?i:i+1][strlen(argv[i])>2?2:0]); p=strtok(mass_range,","); 
	lowest_mass=atoi(p); p=strtok('\0',","); highest_mass=atoi(p);
      }
    if( (argv[i][0]=='-') && (argv[i][1]=='p') ) 
      {
	strcpy(pixels,&argv[strlen(argv[i])>2?i:i+1][strlen(argv[i])>2?2:0]); p=strtok(pixels,","); 
	horizontal_pixels=atoi(p); p=strtok('\0',","); vertical_pixels=atoi(p);
      }
    if( (argv[i][0]=='-') && (argv[i][1]=='B') ) background=atof(&argv[strlen(argv[i])>2?i:i+1][strlen(argv[i])>2?2:0]);
    if( (argv[i][0]=='-') && (argv[i][1]=='S') ) saturation=atof(&argv[strlen(argv[i])>2?i:i+1][strlen(argv[i])>2?2:0]);
    if( strcmp(argv[i],"-noaxes")==0 ) noaxes=1;
  }

  XWIDTH=horizontal_pixels+(1-noaxes)*40; YWIDTH=vertical_pixels+(1-noaxes)*20; /* add pixels for axes */

  printf("checking mzXML dataset %s...",mzXML_filename); fflush(stdout);
  mzXML_file=rampOpenFile(mzXML_filename);
  indexOffset = getIndexOffset (mzXML_file); /* read the index offset */
  pScanIndex = readIndex(mzXML_file, indexOffset, &iLastScan); /* read the scan index into a vector and get LastScan */
  readRunHeader(mzXML_file, pScanIndex, &runHeader, iLastScan);
  printf("done\n"); fflush(stdout);
  
  printf("allocating memory..."); fflush(stdout);
  
  n_scans=last_scan-first_scan+1;
  z=malloc((n_scans>vertical_pixels?n_scans:vertical_pixels)*sizeof(int *));
  if(z==NULL) {fprintf(stderr,"out of memory (terminating)\n"); return -1;}
  stl_height=malloc(vertical_pixels*sizeof(int *));
  if(stl_height==NULL) {fprintf(stderr,"out of memory (terminating)\n"); return -1;}
  smooth=malloc(vertical_pixels*sizeof(int *));
  if(smooth==NULL) {fprintf(stderr,"out of memory (terminating)\n"); return -1;}
  smooth2=malloc(vertical_pixels*sizeof(int *));
  if(smooth2==NULL) {fprintf(stderr,"out of memory (terminating)\n"); return -1;}
  new_spectrum=malloc(horizontal_pixels*sizeof(int));
  scan_number=malloc(n_scans*sizeof(int));
  for(y=0;y<(n_scans>vertical_pixels?n_scans:vertical_pixels);y++) 
    {
      z[y]=malloc(horizontal_pixels*sizeof(int));
      if(z[y]==NULL) {fprintf(stderr,"out of memory (terminating)\n"); return -1;}
    }
 for(y=0;y<vertical_pixels;y++) 
    {
      stl_height[y]=malloc(horizontal_pixels*sizeof(int));
      if(stl_height[y]==NULL) {fprintf(stderr,"out of memory (terminating)\n"); return -1;}
    }
 for(y=0;y<vertical_pixels;y++) 
    {
      smooth[y]=malloc(horizontal_pixels*sizeof(int));
      if(smooth[y]==NULL) {fprintf(stderr,"out of memory (terminating)\n"); return -1;}
    }
 for(y=0;y<vertical_pixels;y++) 
   {
     smooth2[y]=malloc(horizontal_pixels*sizeof(int));
     if(smooth2[y]==NULL) {fprintf(stderr,"out of memory (terminating)\n"); return -1;}
   }
  

  printf("done\nscanning mzXML data for yval range..."); fflush(stdout);

  /* get min_yval and max_yval */

  min_yval=1E10; max_yval=0; y=-1;
  for(scan=first_scan;scan<last_scan;scan++) 
    {
      range[0]=scan;
      readHeader(mzXML_file, pScanIndex[range[0]], &scanHeader);
      if(scanHeader.msLevel==1) /* MS spectrum */
	{
	  y++;
	  scan_number[y]=scan; /* store scan numbers for axis */
	  range[1]=range[0];
	  if(getMsSpect(&mzXML_spectrum,mzXML_file,range)<0) continue; /* MS spectrum not found in mzXML *** */
	  for(i=0;i<mzXML_spectrum.size;i++) /* # of peaks in peaklist in line spectra */
	    {
	      if(mzXML_spectrum.yval[i]>min_yval) min_yval=mzXML_spectrum.yval[i];
	      if(mzXML_spectrum.yval[i]>max_yval) max_yval=mzXML_spectrum.yval[i];
	    }
	}
    }

  printf("done\nparsing mzXML data..."); fflush(stdout);

  if(min_yval>0) log_min_yval=log(min_yval); log_max_yval=log(max_yval);

  /* read in spectrum by spectrum, store in temporary array */
      
  mz_per_pixel=(float)(highest_mass-lowest_mass)/(float)horizontal_pixels; 

  for(y=0;y<vertical_pixels;y++) for(x=0;x<horizontal_pixels;x++) stl_height[y][x]=0;  /* initialize */
  for(y=0;y<(n_scans>vertical_pixels?n_scans:vertical_pixels);y++) for(x=0;x<(highest_mass-lowest_mass)/mz_per_pixel;x++) z[y][x]=0;
    
  y=-1; 

  /* linear - add log and nth root later */
     for(scan=first_scan;scan<last_scan;scan++) 
	{
	  range[0]=scan;
	  readHeader(mzXML_file, pScanIndex[range[0]], &scanHeader);
	  if(scanHeader.msLevel==1) /* MS spectrum */
	    {
	      y++;
	      range[1]=range[0];
	      if(getMsSpect(&mzXML_spectrum,mzXML_file,range)<0) continue; /* MS spectrum not found in mzXML *** */
	      for(i=0;i<mzXML_spectrum.size;i++) /* # of peaks in peaklist in line spectra */
		{
		  if(mzXML_spectrum.xval[i]<lowest_mass) continue;
		  if(mzXML_spectrum.xval[i]>highest_mass) break;
		  if(mzXML_spectrum.yval[i]>background) c=round(255*(mzXML_spectrum.yval[i]-background)/(max_yval*saturation));
		  else c=0;
		  if(c>255) c=255; /* saturated */
		  x=(int)round((mzXML_spectrum.xval[i]-lowest_mass)/mz_per_pixel);
		  if(c>z[y][x]) z[y][x]=c; /* store maximum in each pixel interval */
		}
	    }
	}
 
 
  printf("done\n"); fflush(stdout);

  MS_scans=y; scans_per_pixel=(float)MS_scans/(float)vertical_pixels;



   
  printf("generating STL file from mzXML..."); fflush(stdout);
 
  
  /* plot data */

  if(scans_per_pixel>=1) /* when looking at a whole scene, e.g. LC-MS run with 2000 MS scans = 20 cm with 0.1 mm pixels */
    {
      for(y=0;y<MS_scans;y++)
	{
	  for(x=0;x<horizontal_pixels;x++) stl_height[x+(1-noaxes)*40][YWIDTH-(1-noaxes)*20-(int)round((float)y/(float)scans_per_pixel)]=z[y][x];
	}
    }
  else /* for STL, one would typically expect more than one "pixel" (0.1 mm) per scan for zooms into one isotopic distribution - interpolate here or later? */
    {
      for(y=0;y<MS_scans;y++)
	{
	  for(i=0;i<=ceil(1/scans_per_pixel);i++)
	    {
	      for(x=0;x<horizontal_pixels;x++) stl_height[x+(1-noaxes)*40][YWIDTH-(1-noaxes)*20-i-(int)round((float)y/(float)scans_per_pixel)]=z[y][x];
	    }
	}
    }

  printf("done\nwriting STL...",stl_filename); fflush(stdout); /* add axes later, if at all */

  //xscale=0.7;
  //yscale=0.1;
  //zscale=0.13;
  xscale=1.0;
  yscale=0.12;
  zscale=0.17;

  /*
  if ((outp=fopen(stl_filename, "wb"))==NULL) {printf("\nerror opening output file %s\n",stl_filename); return -1;}
    
  fprintf(outp,"solid ms_data\n");
  for(y=0;y<vertical_pixels-1;y++)
    {
      for(x=1;x<horizontal_pixels;x++)
	{
	  fprintf(outp,"  facet normal 0.0 0.0 0.0\n");
	  fprintf(outp,"    outer loop\n");
	  fprintf(outp,"      vertex %.1f %.1f %.1f\n",(float)x*xscale,(float)y*yscale,(float)stl_height[x][y]*zscale); 
	  fprintf(outp,"      vertex %.1f %.1f %.1f\n",(float)(x-1)*xscale,(float)(y+1)*yscale,(float)stl_height[x-1][y+1]*zscale); 
	  fprintf(outp,"      vertex %.1f %.1f %.1f\n",(float)(x-1)*xscale,(float)y*yscale,(float)stl_height[x-1][y]*zscale);
	  fprintf(outp,"    endloop\n");
	  fprintf(outp,"  endfacet\n");

	  fprintf(outp,"  facet normal 0.0 0.0 0.0\n");
	  fprintf(outp,"    outer loop\n");
	  fprintf(outp,"      vertex %.1f %.1f %.1f\n",(float)x*xscale,(float)y*yscale,(float)stl_height[x][y]*zscale); 
	  fprintf(outp,"      vertex %.1f %.1f %.1f\n",(float)x*xscale,(float)(y+1)*yscale,(float)stl_height[x][y+1]*zscale); 
	  fprintf(outp,"      vertex %.1f %.1f %.1f\n",(float)(x-1)*xscale,(float)(y+1)*yscale,(float)stl_height[x-1][y+1]*zscale);
	  fprintf(outp,"    endloop\n");
	  fprintf(outp,"  endfacet\n");

	} 
    }
  fprintf(outp,"endsolid\n");

  fclose(outp);
 */

  strcpy(stl_filename,strcat(stl_filename,"_smoothed.stl"));
  if ((outp=fopen(stl_filename, "wb"))==NULL) {printf("\nerror opening output file %s\n",stl_filename); return -1;}

  zbase=5.0; /* 5 mm base */
  
  fprintf(outp,"solid ms_data_SG_smoothed\n");

  /* draw box, first bottom, than 4 sides */
  /*
  fprintf(outp,"  facet normal 0.0 0.0 0.0\n"); //bottom 
  fprintf(outp,"    outer loop\n");
  fprintf(outp,"      vertex %.3f %.3f %.3f\n",0.0,0.0,0.0); 
  fprintf(outp,"      vertex %.3f %.3f %.3f\n",0.0,(float)vertical_pixels*yscale,0.0); 
  fprintf(outp,"      vertex %.3f %.3f %.3f\n",(float)horizontal_pixels*xscale,0.0,0.0);
  fprintf(outp,"    endloop\n");
  fprintf(outp,"  endfacet\n");
  fprintf(outp,"  facet normal 0.0 0.0 0.0\n");
  fprintf(outp,"    outer loop\n");
  fprintf(outp,"      vertex %.3f %.3f %.3f\n",0.0,(float)vertical_pixels*yscale,0.0); 
  fprintf(outp,"      vertex %.3f %.3f %.3f\n",(float)horizontal_pixels*xscale,vertical_pixels*yscale,0.0); 
  fprintf(outp,"      vertex %.3f %.3f %.3f\n",(float)horizontal_pixels*xscale,0.0,0.0);
  fprintf(outp,"    endloop\n");
  fprintf(outp,"  endfacet\n");

  fprintf(outp,"  facet normal 0.0 0.0 0.0\n"); //top
  fprintf(outp,"    outer loop\n");
  fprintf(outp,"      vertex %.3f %.3f %.3f\n",0.0,0.0,5.0); 
  fprintf(outp,"      vertex %.3f %.3f %.3f\n",(float)horizontal_pixels*xscale,0.0,5.0);
  fprintf(outp,"      vertex %.3f %.3f %.3f\n",0.0,(float)vertical_pixels*yscale,5.0); 
  fprintf(outp,"    endloop\n");
  fprintf(outp,"  endfacet\n");
  fprintf(outp,"  facet normal 0.0 0.0 0.0\n");
  fprintf(outp,"    outer loop\n");
  fprintf(outp,"      vertex %.3f %.3f %.3f\n",0.0,(float)vertical_pixels*yscale,5.0); 
  fprintf(outp,"      vertex %.3f %.3f %.3f\n",(float)horizontal_pixels*xscale,0.0,5.0);
  fprintf(outp,"      vertex %.3f %.3f %.3f\n",(float)horizontal_pixels*xscale,vertical_pixels*yscale,5.0); 
  fprintf(outp,"    endloop\n");
  fprintf(outp,"  endfacet\n");
  
  fprintf(outp,"  facet normal 0.0 0.0 0.0\n");
  fprintf(outp,"    outer loop\n");
  fprintf(outp,"      vertex %.3f %.3f %.3f\n",0.0,0.0,0.0); 
  fprintf(outp,"      vertex %.3f %.3f %.3f\n",0.0,0.0,5.0); 
  fprintf(outp,"      vertex %.3f %.3f %.3f\n",0.0,(float)vertical_pixels*yscale,0.0);
  fprintf(outp,"    endloop\n");
  fprintf(outp,"  endfacet\n");
  fprintf(outp,"  facet normal 0.0 0.0 0.0\n");
  fprintf(outp,"    outer loop\n");
  fprintf(outp,"      vertex %.3f %.3f %.3f\n",0.0,0.0,5.0); 
  fprintf(outp,"      vertex %.3f %.3f %.3f\n",0.0,(float)vertical_pixels*yscale,5.0); 
  fprintf(outp,"      vertex %.3f %.3f %.3f\n",0.0,(float)vertical_pixels*yscale,0.0);
  fprintf(outp,"    endloop\n");
  fprintf(outp,"  endfacet\n");

  fprintf(outp,"  facet normal 0.0 0.0 0.0\n");
  fprintf(outp,"    outer loop\n");
  fprintf(outp,"      vertex %.3f %.3f %.3f\n",(float)horizontal_pixels*xscale,0.0,0.0); 
  fprintf(outp,"      vertex %.3f %.3f %.3f\n",(float)horizontal_pixels*xscale,(float)vertical_pixels*yscale,0.0);
  fprintf(outp,"      vertex %.3f %.3f %.3f\n",(float)horizontal_pixels*xscale,0.0,5.0); 
  fprintf(outp,"    endloop\n");
  fprintf(outp,"  endfacet\n");
  fprintf(outp,"  facet normal 0.0 0.0 0.0\n");
  fprintf(outp,"    outer loop\n");
  fprintf(outp,"      vertex %.3f %.3f %.3f\n",(float)horizontal_pixels*xscale,0.0,5.0); 
  fprintf(outp,"      vertex %.3f %.3f %.3f\n",(float)horizontal_pixels*xscale,(float)vertical_pixels*yscale,0.0);
  fprintf(outp,"      vertex %.3f %.3f %.3f\n",(float)horizontal_pixels*xscale,(float)vertical_pixels*yscale,5.0); 
  fprintf(outp,"    endloop\n");
  fprintf(outp,"  endfacet\n");

  fprintf(outp,"  facet normal 0.0 0.0 0.0\n");
  fprintf(outp,"    outer loop\n");
  fprintf(outp,"      vertex %.3f %.3f %.3f\n",0.0,0.0,0.0); 
  fprintf(outp,"      vertex %.3f %.3f %.3f\n",(float)horizontal_pixels*xscale,0.0,0.0);
  fprintf(outp,"      vertex %.3f %.3f %.3f\n",0.0,0.0,5.0); 
  fprintf(outp,"    endloop\n");
  fprintf(outp,"  endfacet\n");
  fprintf(outp,"  facet normal 0.0 0.0 0.0\n");
  fprintf(outp,"    outer loop\n");
  fprintf(outp,"      vertex %.3f %.3f %.3f\n",0.0,0.0,5.0); 
  fprintf(outp,"      vertex %.3f %.3f %.3f\n",(float)horizontal_pixels*xscale,0.0,0.0);
  fprintf(outp,"      vertex %.3f %.3f %.3f\n",(float)horizontal_pixels*xscale,0.0,5.0); 
  fprintf(outp,"    endloop\n");
  fprintf(outp,"  endfacet\n");

  fprintf(outp,"  facet normal 0.0 0.0 0.0\n");
  fprintf(outp,"    outer loop\n");
  fprintf(outp,"      vertex %.3f %.3f %.3f\n",0.0,(float)vertical_pixels*yscale,0.0); 
  fprintf(outp,"      vertex %.3f %.3f %.3f\n",0.0,(float)vertical_pixels*yscale,5.0); 
  fprintf(outp,"      vertex %.3f %.3f %.3f\n",(float)horizontal_pixels*xscale,(float)vertical_pixels*yscale,0.0);
  fprintf(outp,"    endloop\n");
  fprintf(outp,"  endfacet\n");
  fprintf(outp,"  facet normal 0.0 0.0 0.0\n");
  fprintf(outp,"    outer loop\n");
  fprintf(outp,"      vertex %.3f %.3f %.3f\n",0.0,(float)vertical_pixels*yscale,5.0); 
  fprintf(outp,"      vertex %.3f %.3f %.3f\n",(float)horizontal_pixels*xscale,(float)vertical_pixels*yscale,5.0); 
  fprintf(outp,"      vertex %.3f %.3f %.3f\n",(float)horizontal_pixels*xscale,(float)vertical_pixels*yscale,0.0);
  fprintf(outp,"    endloop\n");
  fprintf(outp,"  endfacet\n");
  */
  /* absolute box */

  xmin=530; xmax=650;

  box_height=5.0;
  box_base=3.0;
  
  fprintf(outp,"  facet normal 0.0 0.0 0.0\n"); 
  fprintf(outp,"    outer loop\n");
  fprintf(outp,"      vertex %.3f %.3f %.3f\n",(float)xmin*xscale,0.0,box_base); 
  fprintf(outp,"      vertex %.3f %.3f %.3f\n",(float)xmin*xscale,(float)vertical_pixels*yscale,box_base); 
  fprintf(outp,"      vertex %.3f %.3f %.3f\n",(float)xmax*xscale,0.0,box_base);
  fprintf(outp,"    endloop\n");
  fprintf(outp,"  endfacet\n");
  fprintf(outp,"  facet normal 0.0 0.0 0.0\n");
  fprintf(outp,"    outer loop\n");
  fprintf(outp,"      vertex %.3f %.3f %.3f\n",(float)xmin*xscale,(float)vertical_pixels*yscale,box_base); 
  fprintf(outp,"      vertex %.3f %.3f %.3f\n",(float)xmax*xscale,vertical_pixels*yscale,box_base); 
  fprintf(outp,"      vertex %.3f %.3f %.3f\n",(float)xmax*xscale,0.0,box_base);
  fprintf(outp,"    endloop\n");
  fprintf(outp,"  endfacet\n");

  fprintf(outp,"  facet normal 0.0 0.0 0.0\n");
  fprintf(outp,"    outer loop\n");
  fprintf(outp,"      vertex %.3f %.3f %.3f\n",(float)xmin*xscale,0.0,box_height); 
  fprintf(outp,"      vertex %.3f %.3f %.3f\n",(float)xmax*xscale,0.0,box_height);
  fprintf(outp,"      vertex %.3f %.3f %.3f\n",(float)xmin*xscale,(float)vertical_pixels*yscale,box_height); 
  fprintf(outp,"    endloop\n");
  fprintf(outp,"  endfacet\n");
  fprintf(outp,"  facet normal 0.0 0.0 0.0\n");
  fprintf(outp,"    outer loop\n");
  fprintf(outp,"      vertex %.3f %.3f %.3f\n",(float)xmin*xscale,(float)vertical_pixels*yscale,box_height); 
  fprintf(outp,"      vertex %.3f %.3f %.3f\n",(float)xmax*xscale,0.0,box_height);
  fprintf(outp,"      vertex %.3f %.3f %.3f\n",(float)xmax*xscale,vertical_pixels*yscale,box_height); 
  fprintf(outp,"    endloop\n");
  fprintf(outp,"  endfacet\n");
  
  fprintf(outp,"  facet normal 0.0 0.0 0.0\n");
  fprintf(outp,"    outer loop\n");
  fprintf(outp,"      vertex %.3f %.3f %.3f\n",(float)xmin*xscale,0.0,box_base); 
  fprintf(outp,"      vertex %.3f %.3f %.3f\n",(float)xmin*xscale,0.0,box_height); 
  fprintf(outp,"      vertex %.3f %.3f %.3f\n",(float)xmin*xscale,(float)vertical_pixels*yscale,box_base);
  fprintf(outp,"    endloop\n");
  fprintf(outp,"  endfacet\n");
  fprintf(outp,"  facet normal 0.0 0.0 0.0\n");
  fprintf(outp,"    outer loop\n");
  fprintf(outp,"      vertex %.3f %.3f %.3f\n",(float)xmin*xscale,0.0,box_height); 
  fprintf(outp,"      vertex %.3f %.3f %.3f\n",(float)xmin*xscale,(float)vertical_pixels*yscale,box_height); 
  fprintf(outp,"      vertex %.3f %.3f %.3f\n",(float)xmin*xscale,(float)vertical_pixels*yscale,box_base);
  fprintf(outp,"    endloop\n");
  fprintf(outp,"  endfacet\n");

  fprintf(outp,"  facet normal 0.0 0.0 0.0\n");
  fprintf(outp,"    outer loop\n");
  fprintf(outp,"      vertex %.3f %.3f %.3f\n",(float)xmax*xscale,0.0,box_base); 
  fprintf(outp,"      vertex %.3f %.3f %.3f\n",(float)xmax*xscale,(float)vertical_pixels*yscale,box_base);
  fprintf(outp,"      vertex %.3f %.3f %.3f\n",(float)xmax*xscale,0.0,box_height); 
  fprintf(outp,"    endloop\n");
  fprintf(outp,"  endfacet\n");
  fprintf(outp,"  facet normal 0.0 0.0 0.0\n");
  fprintf(outp,"    outer loop\n");
  fprintf(outp,"      vertex %.3f %.3f %.3f\n",(float)xmax*xscale,0.0,box_height); 
  fprintf(outp,"      vertex %.3f %.3f %.3f\n",(float)xmax*xscale,(float)vertical_pixels*yscale,box_base);
  fprintf(outp,"      vertex %.3f %.3f %.3f\n",(float)xmax*xscale,(float)vertical_pixels*yscale,box_height); 
  fprintf(outp,"    endloop\n");
  fprintf(outp,"  endfacet\n");

  fprintf(outp,"  facet normal 0.0 0.0 0.0\n");
  fprintf(outp,"    outer loop\n");
  fprintf(outp,"      vertex %.3f %.3f %.3f\n",(float)xmin*xscale,0.0,box_base); 
  fprintf(outp,"      vertex %.3f %.3f %.3f\n",(float)xmax*xscale,0.0,box_base);
  fprintf(outp,"      vertex %.3f %.3f %.3f\n",(float)xmin*xscale,0.0,box_height); 
  fprintf(outp,"    endloop\n");
  fprintf(outp,"  endfacet\n");
  fprintf(outp,"  facet normal 0.0 0.0 0.0\n");
  fprintf(outp,"    outer loop\n");
  fprintf(outp,"      vertex %.3f %.3f %.3f\n",(float)xmin*xscale,0.0,box_height); 
  fprintf(outp,"      vertex %.3f %.3f %.3f\n",(float)xmax*xscale,0.0,box_base);
  fprintf(outp,"      vertex %.3f %.3f %.3f\n",(float)xmax*xscale,0.0,box_height); 
  fprintf(outp,"    endloop\n");
  fprintf(outp,"  endfacet\n");

  fprintf(outp,"  facet normal 0.0 0.0 0.0\n");
  fprintf(outp,"    outer loop\n");
  fprintf(outp,"      vertex %.3f %.3f %.3f\n",(float)xmin*xscale,(float)vertical_pixels*yscale,box_base); 
  fprintf(outp,"      vertex %.3f %.3f %.3f\n",(float)xmin*xscale,(float)vertical_pixels*yscale,box_height); 
  fprintf(outp,"      vertex %.3f %.3f %.3f\n",(float)xmax*xscale,(float)vertical_pixels*yscale,box_base);
  fprintf(outp,"    endloop\n");
  fprintf(outp,"  endfacet\n");
  fprintf(outp,"  facet normal 0.0 0.0 0.0\n");
  fprintf(outp,"    outer loop\n");
  fprintf(outp,"      vertex %.3f %.3f %.3f\n",(float)xmin*xscale,(float)vertical_pixels*yscale,box_height); 
  fprintf(outp,"      vertex %.3f %.3f %.3f\n",(float)xmax*xscale,(float)vertical_pixels*yscale,box_height); 
  fprintf(outp,"      vertex %.3f %.3f %.3f\n",(float)xmax*xscale,(float)vertical_pixels*yscale,box_base);
  fprintf(outp,"    endloop\n");
  fprintf(outp,"  endfacet\n");






  /*
  // window smoothing  
  for(y=0;y<vertical_pixels;y++)
    {
      for(x=1;x<horizontal_pixels-1;x++)
	{
	  new_spectrum[x]=0.25*stl_height[x-1][y]+0.5*stl_height[x][y]+0.25*stl_height[x+1][y];
	}
      for(x=1;x<horizontal_pixels-1;x++)
	{
	  stl_height[x][y]=new_spectrum[x];
	}
    }
  */

  /* interpolation in zooms k = spikiness factor? */
  /*
  k=4;
  for(y=0;y<vertical_pixels;y++)
    {
      for(x=1;x<horizontal_pixels-k;x++)
	{
	  if( (stl_height[x+1][y]<stl_height[x][y]) && (stl_height[x-1][y]<stl_height[x][y]) ) // x local maximum - interpolate to next local maximum
	    {
	      for(i=x+k;i<horizontal_pixels-1;i++) {if( (stl_height[i-1][y]<stl_height[i][y]) && (stl_height[i+1][y]<stl_height[i][y]) ) break;} // found next local maximum i
	      for(j=x+1;j<i;j++) {stl_height[j][y]=((float)(j-x)/(float)(i-x))*(float)stl_height[i][y]+(1.0-((float)(j-x)/(float)(i-x)))*(float)stl_height[x][y];} // interpolate
	    }
	  x+=k;
	}
    }
  k=4;
  for(y=0;y<vertical_pixels;y++)
    {
      for(x=1;x<horizontal_pixels-k;x++)
	{
	  if( (stl_height[x+1][y]<stl_height[x][y]) && (stl_height[x-1][y]<stl_height[x][y]) ) // x local maximum - interpolate to next local maximum
	    {
	      for(i=x+k;i<horizontal_pixels-1;i++) {if( (stl_height[i-1][y]<stl_height[i][y]) && (stl_height[i+1][y]<stl_height[i][y]) ) break;} // found next local maximum i
	      for(j=x+1;j<i;j++) {stl_height[j][y]=((float)(j-x)/(float)(i-x))*(float)stl_height[i][y]+(1.0-((float)(j-x)/(float)(i-x)))*(float)stl_height[x][y];} // interpolate
	    }
	  x+=k;
	}
    }
  
  k=10;
 for(y=0;y<vertical_pixels;y++)
    {
      for(x=horizontal_pixels-k;x>k+1;x--) // run the same thing but backwards?
	{
	  if( (stl_height[x+1][y]<stl_height[x][y]) && (stl_height[x-1][y]<stl_height[x][y]) ) // x local maximum - interpolate to next local maximum
	    {
	      for(i=x-k;i<horizontal_pixels-1;i++) {if( (stl_height[i-1][y]<stl_height[i][y]) && (stl_height[i+1][y]<stl_height[i][y]) ) break;} // found next local maximum i
	      for(j=x-1;j<i;j++) {stl_height[j][y]=((float)(j-x)/(float)(i-x))*(float)stl_height[i][y]+(1.0-((float)(j-x)/(float)(i-x)))*(float)stl_height[x][y];} // interpolate
	    }
	  x-=k;
	}
    }
  
*/  
  
  
  
  for(y=6;y<vertical_pixels-6;y++)
    {
      for(x=6;x<horizontal_pixels-6;x++) {
	smooth[x][y]=(stl_height[x-6][y]*(-11)+stl_height[x-5][y]*0+stl_height[x-4][y]*9+stl_height[x-3][y]*16+stl_height[x-2][y]*21+stl_height[x-1][y]*24+stl_height[x][y]*25 +    
		      stl_height[x+6][y]*(-11)+stl_height[x+5][y]*0+stl_height[x+4][y]*9+stl_height[x+3][y]*16+stl_height[x+2][y]*21+stl_height[x+1][y]*24)/143;
      }
    }


  /* for(y=5;y<vertical_pixels-5;y++)
    {
      for(x=5;x<horizontal_pixels-5;x++) {
	smooth[x][y]=(stl_height[x-5][y]*(-36)+stl_height[x-4][y]*9+stl_height[x-3][y]*44+stl_height[x-2][y]*69+stl_height[x-1][y]*84+stl_height[x][y]*89 +    
		      stl_height[x+5][y]*(-36)+stl_height[x+4][y]*9+stl_height[x+3][y]*44+stl_height[x+2][y]*69+stl_height[x+1][y]*84)/429;
      }
      }*/

  /*for(y=4;y<vertical_pixels-4;y++)
    {
      for(x=4;x<horizontal_pixels-4;x++) {
	smooth[x][y]=(stl_height[x-4][y]*(-21)+stl_height[x-3][y]*14+stl_height[x-2][y]*39+stl_height[x-1][y]*54+stl_height[x][y]*59 +    
		      stl_height[x+4][y]*(-21)+stl_height[x+3][y]*14+stl_height[x+2][y]*39+stl_height[x+1][y]*54)/231;
      }
      }*/
  
  
  for(x=2;x<horizontal_pixels-2;x++)
    {
      for(y=2;y<vertical_pixels-2;y++) {
	smooth2[x][y]=(smooth[x-2][y-2]*SG[0][0]+smooth[x-2][y-1]*SG[0][1]+smooth[x-2][y]*SG[0][2]+smooth[x-2][y+1]*SG[0][3]+smooth[x-2][y+2]*SG[0][4] +
	               smooth[x-1][y-2]*SG[1][0]+smooth[x-1][y-1]*SG[1][1]+smooth[x-1][y]*SG[1][2]+smooth[x-1][y+1]*SG[1][3]+smooth[x-1][y+2]*SG[1][4] +
	               smooth[x  ][y-2]*SG[2][0]+smooth[x  ][y-1]*SG[2][1]+smooth[x  ][y]*SG[2][2]+smooth[x  ][y+1]*SG[2][3]+smooth[x  ][y+2]*SG[2][4] +
		       smooth[x+1][y-2]*SG[3][0]+smooth[x+1][y-1]*SG[3][1]+smooth[x+1][y]*SG[3][2]+smooth[x+1][y+1]*SG[3][3]+smooth[x+1][y+2]*SG[3][4] +
		       smooth[x+2][y-2]*SG[4][0]+smooth[x+2][y-1]*SG[4][1]+smooth[x+2][y]*SG[4][2]+smooth[x+2][y+1]*SG[4][3]+smooth[x+2][y+2]*SG[4][4]);
      }
    }
  
  
  for(y=0;y<vertical_pixels-1;y++)
    {
      for(x=2;x<horizontal_pixels-1;x++)
	{
	  if( (x>xmin) && (x<=xmax) )
	    {
	      fprintf(outp,"  facet normal 0.0 0.0 0.0\n"); /* in most software, <0,0,0> makes the program calculate the normal using the right hand rule */
	      fprintf(outp,"    outer loop\n");
	      fprintf(outp,"      vertex %.3f %.3f %.3f\n",(float)x*xscale,(float)y*yscale,(float)smooth2[horizontal_pixels-x][y]*zscale+zbase); 
	      fprintf(outp,"      vertex %.3f %.3f %.3f\n",(float)(x-1)*xscale,(float)(y+1)*yscale,(float)smooth2[horizontal_pixels-(x-1)][y+1]*zscale+zbase); 
	      fprintf(outp,"      vertex %.3f %.3f %.3f\n",(float)(x-1)*xscale,(float)y*yscale,(float)smooth2[horizontal_pixels-(x-1)][y]*zscale+zbase);
	      fprintf(outp,"    endloop\n");
	      fprintf(outp,"  endfacet\n");
	      
	      fprintf(outp,"  facet normal 0.0 0.0 0.0\n"); /* in most software, <0,0,0> makes the program calculate the normal using the right hand rule */
	      fprintf(outp,"    outer loop\n");
	      fprintf(outp,"      vertex %.3f %.3f %.3f\n",(float)x*xscale,(float)y*yscale,(float)smooth2[horizontal_pixels-x][y]*zscale+zbase); 
	      fprintf(outp,"      vertex %.3f %.3f %.3f\n",(float)x*xscale,(float)(y+1)*yscale,(float)smooth2[horizontal_pixels-x][y+1]*zscale+zbase); 
	      fprintf(outp,"      vertex %.3f %.3f %.3f\n",(float)(x-1)*xscale,(float)(y+1)*yscale,(float)smooth2[horizontal_pixels-(x-1)][y+1]*zscale+zbase);
	      fprintf(outp,"    endloop\n");
	      fprintf(outp,"  endfacet\n");
	    }

	} 
    }
  fprintf(outp,"endsolid\n");

  fclose(outp);

  printf("done (STL saved to '%s')\n",stl_filename); fflush(stdout);
  
  return 0;
}

