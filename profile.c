/* 
** profile.c
**
** Program written in order to calculate profile
**
** written by Marcel Zemp
*/

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <malloc.h>
#include <assert.h>
#include <iof.h>

typedef struct profilestructure {

    int Ntot;
    int Ngas;
    int Ndark;
    int Nstar;
    int Nenctot;
    int Nencgas;
    int Nencdark;
    int Nencstar;
    double ri;
    double rm;
    double ro;
    double vol;
    double Mtot;
    double Mgas;
    double Mdark;
    double Mstar;
    double Menctot;
    double Mencgas;
    double Mencdark;
    double Mencstar;
    double veltot[3];
    double vel2tot[6];
    double velgas[3];
    double vel2gas[6];
    double veldark[3];
    double vel2dark[6];
    double velstar[3];
    double vel2star[6];
    double Ltot[3];
    double Lgas[3];
    double Ldark[3];
    double Lstar[3];
    } PS;

typedef struct halodata {

    int ID;
    double rcentre[3];
    double vcentre[3];
    double rbg, Mbg;
    double rcrit, Mcrit;
    double rvcmax, Mrvcmax, vcmax;
    PS *ps;
    } HD;

typedef struct generalinfo {

    int positionprecision;
    int gridtype;
    int projectionvariant;
    int centretype;
    int readhalocataloguefile;
    int rmaxfromhalocatalogue;
    int dataformat;
    int halocatalogueformat;
    int verboselevel;
    int doswap;
    int Nbin;
    int NHalo, Ngas, Ndark, Nstar;
    int NParticlesPerBlock, NParticlesInBlock, NBlock;
    double acurrent;
    double rhoencbg, rhoenccrit;
    double rhocrit0, OmegaM0, OmegaK0, OmegaL0, OmegaM, rhocrit, E, Deltabg, Deltacrit;
    double rmin, rmax;
    double G;
    double bc[6];
    double binfactor;
    char outputname[256];
    } GI;

void calculate_unit_vectors(double pos[3], double erad[3], double ephi[3], double etheta[3]) {

    double dist;
    double cosphi, sinphi, costheta, sintheta;

    /*
    ** Calculate cosphi & sinphi
    */
    dist = sqrt(pos[0]*pos[0]+pos[1]*pos[1]);
    cosphi = pos[0]/dist;
    sinphi = pos[1]/dist;
    if ((pos[0] == 0) && (pos[1] == 0)) {
        cosphi = 1;
        sinphi = 0;
        }
    /*
    ** Calculate costheta & sintheta
    */
    dist = sqrt(pos[0]*pos[0]+pos[1]*pos[1]+pos[2]*pos[2]);
    costheta = pos[2]/dist;
    sintheta = sqrt(1-costheta*costheta);
    /*
    ** Calculate unit vectors
    */
    erad[0] = sintheta*cosphi;
    erad[1] = sintheta*sinphi;
    erad[2] = costheta;
    ephi[0] = -sinphi;
    ephi[1] = cosphi;
    ephi[2] = 0;
    etheta[0] = -costheta*cosphi;
    etheta[1] = -costheta*sinphi;
    etheta[2] = sintheta;
    }

double correct_position(double c, double r, double l) {

    double d;

    d = r-c;
    if (d > 0.5*l) return r - l;
    else if (d < -0.5*l) return r + l;
    else return r;
    }

double put_in_box(double r, double lmin, double lmax) {

    double l;

    l = lmax - lmin;
    if (r < lmin) return r + l;
    else if (r > lmax) return r - l;
    else return r;
    }

double Ecosmo(double a, double OmegaM0, double OmegaL0, double OmegaK0) {

    return sqrt(OmegaM0*pow(a,-3.0) + OmegaL0 + OmegaK0*pow(a,-2.0));
    }

void read_art_header(FILE *fp, char *ab, ART_HEADER *ah, int *doswap) {

    int header, trailer;

    *doswap = 0;
    assert(fread(&header,sizeof(int),1,fp) == 1);
    if (header != 45+sizeof(ART_HEADER)) {
	*doswap = 1;
	flip_4byte(&header,sizeof(int),1);
	}
    assert(header == 45+sizeof(ART_HEADER));
    assert(fread(ab,sizeof(char),45,fp) == 45);
    assert(fread(ah,sizeof(ART_HEADER),1,fp) == 1);
    if (*doswap) flip_4byte(ah,sizeof(ART_HEADER),1);
    assert(fread(&trailer,sizeof(int),1,fp) == 1);
    if (*doswap) flip_4byte(&trailer,sizeof(int),1);
    assert(header == trailer);
    }

void usage(void);
void read_halocatalogue_6DFOF(FILE (*), int, int, int, double, double, int, double, double, double (*), HD (*), int (*));
void put_dp_in_bins(HD (*), DARK_PARTICLE (*), GI);
void calculate_halo_properties(HD (*), GI);
void write_output(HD (*), GI);

int main(int argc, char **argv) {

    int i, j;
    char halocataloguefilename[256];
    char HeaderFileName[256], DataFileName[256], banner[45];
    GI gi;
    TIPSY_HEADER th;
    GAS_PARTICLE *gp = NULL;
    DARK_PARTICLE *dp = NULL;
    STAR_PARTICLE *sp = NULL;
    GAS_PARTICLE_DPP *gpdpp = NULL;
    DARK_PARTICLE_DPP *dpdpp = NULL;
    STAR_PARTICLE_DPP *spdpp = NULL;
    ART_HEADER ah;
    HD *hd = NULL;
    FILE *halocataloguefile;
    XDR xdrs;
    FILE *HeaderFile;
    FILE *PosXFile, *PosYFile, *PosZFile;
    FILE *VelXFile, *VelYFile, *VelZFile;

    /*
    ** Set some default values
    */

    gi.rhocrit0 = 1;
    gi.OmegaM0 = 0.3;
    gi.OmegaL0 = 0.7;
    gi.OmegaK0 = 0;
    gi.Deltabg = 200;
    gi.Deltacrit = 0;
    gi.G = 1;
    gi.rmin = 0;
    gi.rmax = 1;
    gi.Nbin = 0;
    gi.binfactor = 5;

    /*
    ** Read in arguments
    */

    i = 1;
    while (i < argc) {
	if (strcmp(argv[i],"-spp") == 0) {
            gi.positionprecision = 0;
            i++;
            }
        else if (strcmp(argv[i],"-dpp") == 0) {
            gi.positionprecision = 1;
            i++;
            }
        else if (strcmp(argv[i],"-lin") == 0) {
	    gi.gridtype = 0;
            i++;
            }
        else if (strcmp(argv[i],"-log") == 0) {
	    gi.gridtype = 1;
            i++;
            }
        else if (strcmp(argv[i],"-v") == 0) {
	    gi.verboselevel = 1;
            i++;
            }
        else if (strcmp(argv[i],"-pv") == 0) {
            i++;
            if (i >= argc) usage();
            gi.projectionvariant = atoi(argv[i]);
            i++;
            }
	else if (strcmp(argv[i],"-rmin") == 0) {
	    i++;
            if (i >= argc) usage();
	    gi.rmin = atof(argv[i]);
	    i++;
	    }
	else if (strcmp(argv[i],"-rmax") == 0) {
	    i++;
            if (i >= argc) usage();
	    gi.rmax = atof(argv[i]);
	    i++;
	    }
	else if (strcmp(argv[i],"-binfactor") == 0) {
	    i++;
            if (i >= argc) usage();
	    gi.binfactor = atof(argv[i]);
	    i++;
	    }
	else if (strcmp(argv[i],"-Nbin") == 0) {
	    i++;
            if (i >= argc) usage();
	    gi.Nbin = atof(argv[i]);
	    i++;
	    }
        else if (strcmp(argv[i],"-Delta_bg") == 0) {
            i++;
            if (i >= argc) usage();
            gi.Deltabg = atof(argv[i]);
            i++;
            }
        else if (strcmp(argv[i],"-Delta_crit") == 0) {
            i++;
            if (i >= argc) usage();
            gi.Deltacrit = atof(argv[i]);
            i++;
            }
        else if (strcmp(argv[i],"-OmegaM0") == 0) {
            i++;
            if (i >= argc) usage();
            gi.OmegaM0 = atof(argv[i]);
            i++;
            }
        else if (strcmp(argv[i],"-OmegaK0") == 0) {
            i++;
            if (i >= argc) usage();
            gi.OmegaK0 = atof(argv[i]);
            i++;
            }
        else if (strcmp(argv[i],"-OmegaL0") == 0) {
            i++;
            if (i >= argc) usage();
            gi.OmegaL0 = atof(argv[i]);
            i++;
            }
        else if (strcmp(argv[i],"-rhocrit0") == 0) {
            i++;
            if (i >= argc) usage();
            gi.rhocrit0 = atof(argv[i]);
            i++;
            }
        else if (strcmp(argv[i],"-G") == 0) {
            i++;
            if (i >= argc) usage();
            gi.G = atof(argv[i]);
            i++;
            }
        else if (strcmp(argv[i],"-com") == 0) {
            gi.centretype = 0;
            i++;
            }
        else if (strcmp(argv[i],"-potmin") == 0) {
            gi.centretype = 1;
            i++;
            }
        else if (strcmp(argv[i],"-denmax") == 0) {
            gi.centretype = 2;
            i++;
            }
        else if (strcmp(argv[i],"-rmaxfromhalocatalogue") == 0) {
            gi.rmaxfromhalocatalogue = 1;
            i++;
            }
        else if (strcmp(argv[i],"-dataformat") == 0) {
            i++;
            if (i >= argc) usage();
	    gi.dataformat = atoi(argv[i]);
            i++;
            }
        else if (strcmp(argv[i],"-halocatalogueformat") == 0) {
            i++;
            if (i >= argc) usage();
	    gi.halocatalogueformat = atoi(argv[i]);
            i++;
            }
        else if (strcmp(argv[i],"-halocatalogue") == 0) {
	    gi.readhalocataloguefile = 1;
            i++;
            if (i >= argc) usage();
            strcpy(halocataloguefilename,argv[i]);
            i++;
            }
        else if (strcmp(argv[i],"-output") == 0) {
            i++;
            if (i >= argc) usage();
            strcpy(gi.outputname,argv[i]);
            i++;
            }
        else if (strcmp(argv[i],"-header") == 0) {
            i++;
            if (i >= argc) usage();
            strcpy(HeaderFileName,argv[i]);
            i++;
            }
        else if (strcmp(argv[i],"-data") == 0) {
            i++;
            if (i >= argc) usage();
            strcpy(DataFileName,argv[i]);
            i++;
            }
	else if ((strcmp(argv[i],"-h") == 0) || (strcmp(argv[i],"-help") == 0)) {
	    usage();
	    }
	else {
	    usage();
	    }
	}

    /*
    ** Read header files
    */

    if (gi.dataformat == 0) {
	xdrstdio_create(&xdrs,stdin,XDR_DECODE);
	read_tipsy_standard_header(&xdrs,&th);
	gi.acurrent = th.time;
	gi.Ngas = th.ngas;
	gi.Ndark = th.ndark;
	gi.Nstar = th.nstar;
	gi.NParticlesPerBlock = 1000000;
	gi.bc[0] = -0.5;
	gi.bc[1] = -0.5;
	gi.bc[2] = -0.5;
	gi.bc[3] = 0.5;
	gi.bc[4] = 0.5;
	gi.bc[5] = 0.5;
	}
    else if (gi.dataformat == 1) {
	HeaderFile = fopen(HeaderFileName,"r");
	assert(HeaderFile != NULL);
	read_art_header(HeaderFile,banner,&ah,&gi.doswap);
	fclose(HeaderFile);
	gi.acurrent = ah.aunin;
	gi.OmegaM0 = ah.OmM0;
	gi.OmegaL0 = ah.OmL0;
	gi.OmegaK0 = ah.OmK0;
	gi.Ngas = 0;
	gi.Ndark = 0;
	gi.Nstar = 0;
	gi.NParticlesPerBlock = ah.Nrow*ah.Nrow;
	gi.bc[0] = 0;
	gi.bc[1] = 0;
	gi.bc[2] = 0;
	gi.bc[3] = 1;
	gi.bc[4] = 1;
	gi.bc[5] = 1;
	}
    else {
	fprintf(stderr,"Not supported format!\n");
	exit(1);
	}

    fprintf(stderr,"banner: %s\n",banner);
    fprintf(stderr,"a = %g doswap %d \n",gi.acurrent,gi.doswap);
//    exit(1);

    /*
    ** Calculate cosmology relevant stuff
    ** Densities in comoving coordinates 
    */

    gi.rhoencbg = gi.Deltabg*gi.OmegaM0*gi.rhocrit0;
    gi.E = Ecosmo(gi.acurrent,gi.OmegaM0,gi.OmegaL0,gi.OmegaK0);
    gi.OmegaM = gi.OmegaM0/(pow(gi.acurrent,3)*gi.E*gi.E);
    gi.rhocrit = gi.rhocrit0*gi.E*gi.E;
    if (gi.Deltacrit == 0) {
	if (gi.OmegaK0 == 0) {
	    gi.Deltacrit = 178*pow(gi.OmegaM,0.45);
	    }
	else if (gi.OmegaL0 == 0) {
	    gi.Deltacrit = 178*pow(gi.OmegaM,0.3);
	    }
	}
    gi.rhoenccrit = gi.Deltacrit*gi.rhocrit*pow(gi.acurrent,3);

    /*
    ** Read halo catalogue
    */

    // initialise HD
    halocataloguefile = fopen(halocataloguefilename,"r");
    assert(halocataloguefile != NULL);
    if (gi.halocatalogueformat == 0) {
	read_halocatalogue_6DFOF(halocataloguefile,gi.centretype,gi.gridtype,gi.rmaxfromhalocatalogue,gi.rmin,gi.rmax,gi.Nbin,gi.rhoencbg,gi.binfactor,gi.bc,hd,&gi.NHalo);
	}
    else if (gi.halocatalogueformat == 1) {
//	read_halocatalogue_HFIND();
	}
    fclose(halocataloguefile);

    /*
    ** Harvest data
    */
    if (gi.dataformat == 0) {
	if (gi.positionprecision == 0) {
	    /*
	    ** Dark Matter
	    */
	    dp = realloc(dp,gi.NParticlesPerBlock*sizeof(DARK_PARTICLE));
	    gi.NBlock = (gi.Ndark+gi.NParticlesPerBlock-1)/gi.NParticlesPerBlock;
	    for (i = 0; i < gi.NBlock; i++) {
		if (i == gi.NBlock-1) {
		    gi.NParticlesInBlock = gi.Ndark-((gi.NBlock-1)*gi.NParticlesPerBlock);
		    }
		else {
		    gi.NParticlesInBlock = gi.NParticlesPerBlock;
		    }
		for (j = 0; j < gi.NParticlesInBlock; j++) {
		    read_tipsy_standard_dark(&xdrs,&dp[j]);
		    }
		put_dp_in_bins(hd,dp,gi);
		}
	    free(dp);
	    }
	}
    else if (gi.dataformat == 1) {
	
	}

    /*
    ** Calculate halo properties
    */
 
    calculate_halo_properties(hd,gi); 

    /*
    ** Write output
    */

    write_output(hd,gi);

    /*
    ** Some more output if desired
    */

/*
    if (gi.verboselevel >= 1) {
        fprintf(stderr,"Used values:\n");
        fprintf(stderr,"Delta_bg    : %.6e\n",Deltabg);
        fprintf(stderr,"Delta_crit  : %.6e\n",Deltacrit);
	fprintf(stderr,"rhoenc_bg   : %.6e MU LU^{-3} (comoving)\n",rhoencbg);
	fprintf(stderr,"rhoenc_crit : %.6e MU LU^{-3} (comoving)\n",rhoenccrit);
        fprintf(stderr,"rhocrit0    : %.6e MU LU^{-3}\n",rhocrit0);
        fprintf(stderr,"OmegaM0     : %.6e\n",OmegaM0);
        fprintf(stderr,"OmegaL0     : %.6e\n",OmegaL0);
        fprintf(stderr,"OmegaK0     : %.6e\n",OmegaK0);
        fprintf(stderr,"G           : %.6e LU^3 MU^{-1} TU^{-1}\n",G);
	fprintf(stderr,"a           : %.6e\n",acurrent);
	fprintf(stderr,"E           : %.6e\n",E);
        fprintf(stderr,"binfactor   : %.6e\n",binfactor);
        }
*/
    exit(0);

    }

void usage(void) {

    fprintf(stderr,"\n");
    fprintf(stderr,"Program calculates the profile of the input file around r_cen.\n");
    fprintf(stderr,"\n");
    fprintf(stderr,"You can specify the following arguments:\n");
    fprintf(stderr,"\n");
    fprintf(stderr,"-spp                : set this flag if input and output files have single precision positions (default)\n");
    fprintf(stderr,"-dpp                : set this flag if input and output files have double precision positions\n");
    fprintf(stderr,"-lin                : set this flag for linear grid\n");
    fprintf(stderr,"-log                : set this flag for logarithmic grid (default)\n");
    fprintf(stderr,"-pv <value>         : projection variant - 0:along axis, 1:spherical coordinates (default: 0)\n");
    fprintf(stderr,"-rmin <value>       : minimum grid radius [LU] (default 0 LU)\n");
    fprintf(stderr,"-rmax <value>       : maximum grid radius [LU] (default 1 LU)\n");
    fprintf(stderr,"-Nbin <value>       : number of bins between rmin and rmax\n");
    fprintf(stderr,"-N <value>          : number of haloes (default: 1)\n");
    fprintf(stderr,"-rxcen <value>      : x-coordinate of centre [LU] (default: 0 LU)\n");
    fprintf(stderr,"-rycen <value>      : y-coordinate of centre [LU] (default: 0 LU)\n");
    fprintf(stderr,"-rzcen <value>      : z-coordinate of centre [LU] (default: 0 LU)\n");
    fprintf(stderr,"-vxcen <value>      : x-velocity of centre [VU] (default: 0 VU)\n");
    fprintf(stderr,"-vycen <value>      : y-velocity of centre [VU] (default: 0 VU)\n");
    fprintf(stderr,"-vzcen <value>      : z-velocity of centre [VU] (default: 0 VU)\n");
    fprintf(stderr,"-stats <name>       : name of stats file\n");
    fprintf(stderr,"-com                : set this flag for centre-of-mass centres from stats file\n");
    fprintf(stderr,"-potmin             : set this flag for potmin centres from stats file\n");
    fprintf(stderr,"-denmax             : set this flag for denmax centres from stats file (default)\n");
    fprintf(stderr,"-rmaxfromhalocatalogue      : set this flag for rmax determined from stats file\n");
    fprintf(stderr,"-binfactor <value>  : extra factor for rmax determined form stats file (default: 5)\n");
    fprintf(stderr,"-rhocrit0 <value>   : critical density of the universe [MU LU^{-3}] (default: 1 MU LU^{-3})\n");
    fprintf(stderr,"-OmegaM0 <value>    : OmegaM0 value (default: 0.3)\n");
    fprintf(stderr,"-OmegaL0 <value>    : OmegaL0 value (default: 0.7)\n");
    fprintf(stderr,"-OmegaK0 <value>    : OmegaK0 value (default: 0.0)\n");
    fprintf(stderr,"-Delta_bg <value>   : overdensity with respect to background density (default: 200)\n");
    fprintf(stderr,"-Delta_crit <value> : overdensity with respect to critical density (default: 178*(OmegaM^0.45) [OmegaK0=0] / 178*(OmegaM^0.3) [OmegaL0=0])\n");
    fprintf(stderr,"-G <value>          : gravitational constant [LU^3 MU^{-1} TU^{-1}] (default: 1 LU^3 MU^{-1} TU^{-1})\n");
    fprintf(stderr,"-output <name>      : name of output files (endings .statistics and .profiles appended)\n");
    fprintf(stderr,"-v                  : more informative output to screen\n");
    fprintf(stderr,"< <name>            : name of input file in tipsy standard binary format\n");
    fprintf(stderr,"\n");
    exit(1);
    }

void read_halocatalogue_6DFOF(FILE *halocataloguefile, int centretype, int gridtype, int rmaxfromhalocatalogue, double rmin, double rmax, int Nbin, double rhoencbg, double binfactor, double bc[6], HD *hd, int *NHalo) {

    int SizeHaloData = 10000;
    int i, j, k, ID, N, idummy, NHaloRead;
    float fdummy;
    double DarkMass, radius1, radius2, vd1D;
    double rxcom, rycom, rzcom, rxpotmin, rypotmin, rzpotmin, rxdenmax, rydenmax, rzdenmax, vx, vy, vz;
    double rmaxestimate, dr;

    assert(Nbin > 0);
    assert(rmin >= 0);
    assert(rmax >= 0);
    assert(rmax > rmin);
    dr = 0;
    NHaloRead = 0;
    hd = realloc(hd,SizeHaloData*sizeof(HD));
    assert(hd != NULL);
    for (i = 0; i < SizeHaloData; i++) {
	hd[i].ps = realloc(hd[i].ps,(Nbin+1)*sizeof(PS));
	assert(hd[i].ps != NULL);
	}
    while (1) {
	fscanf(halocataloguefile,"%i",&idummy); ID = idummy;
	fscanf(halocataloguefile,"%i",&idummy); N = idummy;
	fscanf(halocataloguefile,"%g",&fdummy); DarkMass = fdummy;
	fscanf(halocataloguefile,"%g",&fdummy);
	fscanf(halocataloguefile,"%g",&fdummy);
	fscanf(halocataloguefile,"%g",&fdummy);
	fscanf(halocataloguefile,"%g",&fdummy); radius1 = fdummy; /* (sum_j rmax[j]-rmin[j])/6 */
	fscanf(halocataloguefile,"%g",&fdummy); radius2 = fdummy; /* dispersion in coordinates */
	fscanf(halocataloguefile,"%g",&fdummy); vd1D = fdummy;
	fscanf(halocataloguefile,"%g",&fdummy);
	fscanf(halocataloguefile,"%g",&fdummy);
	fscanf(halocataloguefile,"%g",&fdummy);
	fscanf(halocataloguefile,"%g",&fdummy); rxcom = put_in_box(fdummy,bc[0],bc[3]);
	fscanf(halocataloguefile,"%g",&fdummy); rycom = put_in_box(fdummy,bc[1],bc[4]);
	fscanf(halocataloguefile,"%g",&fdummy); rzcom = put_in_box(fdummy,bc[2],bc[5]);
	fscanf(halocataloguefile,"%g",&fdummy); rxpotmin = put_in_box(fdummy,bc[0],bc[3]);
	fscanf(halocataloguefile,"%g",&fdummy); rypotmin = put_in_box(fdummy,bc[1],bc[4]);
	fscanf(halocataloguefile,"%g",&fdummy); rzpotmin = put_in_box(fdummy,bc[2],bc[5]);
	fscanf(halocataloguefile,"%g",&fdummy); rxdenmax = put_in_box(fdummy,bc[0],bc[3]);
	fscanf(halocataloguefile,"%g",&fdummy); rydenmax = put_in_box(fdummy,bc[1],bc[4]);
	fscanf(halocataloguefile,"%g",&fdummy); rzdenmax = put_in_box(fdummy,bc[2],bc[5]);
	fscanf(halocataloguefile,"%g",&fdummy); vx = fdummy;
	fscanf(halocataloguefile,"%g",&fdummy); vy = fdummy;
	fscanf(halocataloguefile,"%g",&fdummy); vz = fdummy;
	fscanf(halocataloguefile,"%g",&fdummy);
	fscanf(halocataloguefile,"%g",&fdummy);
	fscanf(halocataloguefile,"%g",&fdummy);
	fscanf(halocataloguefile,"%g",&fdummy);
	fscanf(halocataloguefile,"%g",&fdummy);
	if (feof(halocataloguefile)) break;
	NHaloRead++;
	if (SizeHaloData < NHaloRead){
	    SizeHaloData += 10000; 
	    hd = realloc(hd,SizeHaloData*sizeof(HD));
	    assert(hd != NULL);
	    for (i = 0; i < SizeHaloData; i++) {
		hd[i].ps = realloc(hd[i].ps,(Nbin+1)*sizeof(PS));
		assert(hd[i].ps != NULL);
		}
	    }
	i = NHaloRead-1;
	hd[i].ID = ID;
	if (centretype == 0) {
	    hd[i].rcentre[0] = rxcom;
	    hd[i].rcentre[1] = rycom;
	    hd[i].rcentre[2] = rzcom;
	    }
	else if (centretype == 1) {
	    hd[i].rcentre[0] = rxpotmin;
	    hd[i].rcentre[1] = rypotmin;
	    hd[i].rcentre[2] = rzpotmin;
	    }
	else if (centretype == 2) {
	    hd[i].rcentre[0] = rxdenmax;
	    hd[i].rcentre[1] = rydenmax;
	    hd[i].rcentre[2] = rzdenmax;
	    }
	hd[i].vcentre[0] = vx;
	hd[i].vcentre[1] = vy;
	hd[i].vcentre[2] = vz;
	hd[i].rbg = 0;
	hd[i].Mbg = 0;
	hd[i].rcrit = 0;
	hd[i].Mcrit = 0;
	hd[i].rvcmax = 0;
	hd[i].Mrvcmax = 0;
	hd[i].vcmax = 0;
	if (rmaxfromhalocatalogue == 1) {
	    /*
	    ** Estimate maximum radius; assume isothermal sphere scaling 
	    */
	    rmaxestimate = sqrt((3*DarkMass/(4*M_PI*radius2*radius2*radius2))/rhoencbg)*radius2*binfactor;
	    assert(rmaxestimate > 0);
	    if (gridtype == 0) {
		dr = (rmaxestimate-rmin)/Nbin;
		}
	    else if (gridtype == 1) {
		dr = (log(rmaxestimate)-log(rmin))/Nbin;
		}
	    }
	else {
	    if (gridtype == 0) {
		dr = (rmax-rmin)/Nbin;
		}
	    else if (gridtype == 1) {
		dr = (log(rmax)-log(rmin))/Nbin;
		}
	    }
	for (j = 1; j < (Nbin+1); j++) {
	    if (gridtype == 0) {
		hd[i].ps[j].ri = rmin + (j-1)*dr;
		hd[i].ps[j].ro = rmin + j*dr;
		}
	    else if (gridtype == 1) {
		hd[i].ps[j].ri = exp(log(rmin) + (j-1)*dr);
		hd[i].ps[j].ro = exp(log(rmin) + j*dr);
		}
	    }
	hd[i].ps[j].Ntot = 0;
	hd[i].ps[j].Ngas = 0;
	hd[i].ps[j].Ndark = 0;
	hd[i].ps[j].Nstar = 0;
	hd[i].ps[j].Nenctot = 0;
	hd[i].ps[j].Nencgas = 0;
	hd[i].ps[j].Nencdark = 0;
	hd[i].ps[j].Nencstar = 0;
	hd[i].ps[j].Mtot = 0;
	hd[i].ps[j].Mgas = 0;
	hd[i].ps[j].Mdark = 0;
	hd[i].ps[j].Mstar = 0;
	hd[i].ps[j].Menctot = 0;
	hd[i].ps[j].Mencgas = 0;
	hd[i].ps[j].Mencdark = 0;
	hd[i].ps[j].Mencstar = 0;
	for (k = 0; k < 3; k++) {
	    hd[i].ps[j].veltot[k] = 0;
	    hd[i].ps[j].velgas[k] = 0;
	    hd[i].ps[j].veldark[k] = 0;
	    hd[i].ps[j].velstar[k] = 0;
	    hd[i].ps[j].Ltot[k] = 0;
	    hd[i].ps[j].Lgas[k] = 0;
	    hd[i].ps[j].Ldark[k] = 0;
	    hd[i].ps[j].Lstar[k] = 0;
	    }
	for (k = 0; k < 6; k++) {
	    hd[i].ps[j].vel2tot[k] = 0;
	    hd[i].ps[j].vel2gas[k] = 0;
	    hd[i].ps[j].vel2dark[k] = 0;
	    hd[i].ps[j].vel2star[k] = 0;
	    }
	}
    *NHalo = NHaloRead;
    }

void put_dp_in_bins(HD *hd, DARK_PARTICLE *dp, GI gi) {

    int i, j, k, l;
    double pos[3], vel[3], velproj[3];
    double erad[3], ephi[3], etheta[3];
    double r;

    for (i = 0; i < gi.NParticlesPerBlock; i++) {
	for (j = 0; j < gi.NHalo; j++) {
	    for (k = 0; k < 3; k++) {
		pos[k] = correct_position(hd[j].rcentre[k],dp[i].pos[k],gi.bc[k+3]-gi.bc[k]);
		pos[k] = pos[k]-hd[j].rcentre[k];
		}
	    r = sqrt(pos[0]*pos[0]+pos[1]*pos[1]+pos[2]*pos[2]);
	    if (r <= hd[j].ps[gi.Nbin].ro) {
		for (k = 0; k < 3; k++) {
		    vel[k] = dp[i].vel[k]-hd[j].vcentre[k];
		    }
		/*
		** Go trough bins from outside => larger bin volume further out
		*/
		for (l = gi.Nbin; l >=0; l--) {
		    if ((hd[j].ps[l].ri <= r) && (hd[j].ps[l].ro > r)) {
			if (gi.projectionvariant == 0) {
			    velproj[0] = vel[0];
			    velproj[1] = vel[1];
			    velproj[2] = vel[2];
			    }
			else if (gi.projectionvariant == 1) {
			    calculate_unit_vectors(pos,erad,ephi,etheta);
			    velproj[0] = vel[0]*erad[0]+vel[1]*erad[1]+vel[2]*erad[2];
			    velproj[1] = vel[0]*ephi[0]+vel[1]*ephi[1]+vel[2]*ephi[2];
			    velproj[2] = vel[0]*etheta[0]+vel[1]*etheta[1]+vel[2]*etheta[2];
			    }
			hd[j].ps[l].Ntot++;
			hd[j].ps[l].Ndark++;
			hd[j].ps[l].Mtot += dp[i].mass;
			hd[j].ps[l].Mdark += dp[i].mass;
			for (k = 0; k < 3; k++) {
			    hd[j].ps[l].veltot[k] += velproj[k];
			    hd[j].ps[l].veldark[k] += velproj[k];
			    hd[j].ps[l].vel2tot[k] += velproj[k]*velproj[k];
			    hd[j].ps[l].vel2dark[k] += velproj[k]*velproj[k];
			    }
			hd[j].ps[l].vel2tot[3] += velproj[0]*velproj[1];
			hd[j].ps[l].vel2tot[4] += velproj[0]*velproj[2];
			hd[j].ps[l].vel2tot[5] += velproj[1]*velproj[2];
			hd[j].ps[l].vel2dark[3] += velproj[0]*velproj[1];
			hd[j].ps[l].vel2dark[4] += velproj[0]*velproj[2];
			hd[j].ps[l].vel2dark[5] += velproj[1]*velproj[2];
			hd[j].ps[l].Ltot[0] += dp[i].mass*(pos[1]*vel[2] - pos[2]*vel[1]);
			hd[j].ps[l].Ltot[1] += dp[i].mass*(pos[2]*vel[0] - pos[0]*vel[2]);
			hd[j].ps[l].Ltot[2] += dp[i].mass*(pos[0]*vel[1] - pos[1]*vel[0]);
			hd[j].ps[l].Ldark[0] += dp[i].mass*(pos[1]*vel[2] - pos[2]*vel[1]);
			hd[j].ps[l].Ldark[1] += dp[i].mass*(pos[2]*vel[0] - pos[0]*vel[2]);
			hd[j].ps[l].Ldark[2] += dp[i].mass*(pos[0]*vel[1] - pos[1]*vel[0]);
			break;
			}
		    }
		}
	    }
	}
    }

void calculate_halo_properties(HD *hd, GI gi) {

    int i, j, k;
    double radius[2], rhoenc[2], Menc[2];
    double m, d;
    
    for (i = 0; i < gi.NHalo; i++) {
	
	/*
	** Calculate derived properties
	*/
	
	for (j = 0; j < (gi.Nbin+1); j++) {
	    hd[i].ps[j].vol = 4*M_PI*(hd[i].ps[j].ro*hd[i].ps[j].ro*hd[i].ps[j].ro - hd[i].ps[j].ri*hd[i].ps[j].ri*hd[i].ps[j].ri)/3.0;
	    if (gi.gridtype == 0) {
		hd[i].ps[j].rm = (hd[i].ps[j].ri+hd[i].ps[j].ro)/2.0;
		}
	    else if (gi.gridtype == 1) {
		if (j == 0) {
		    hd[i].ps[j].rm = exp((3*log(hd[i].ps[j+1].ri)-log(hd[i].ps[j+1].ro))/2.0);
		    }
		else {
		    hd[i].ps[j].rm = exp((log(hd[i].ps[j].ri)+log(hd[i].ps[j].ro))/2.0);
		    }
		}
	    for (k = 0; k <= j; k++) {
		hd[i].ps[j].Nenctot  += hd[i].ps[k].Ntot;
		hd[i].ps[j].Nencgas  += hd[i].ps[k].Ngas;
		hd[i].ps[j].Nencdark += hd[i].ps[k].Ndark;
		hd[i].ps[j].Nencstar += hd[i].ps[k].Nstar;
		hd[i].ps[j].Menctot  += hd[i].ps[k].Mtot;
		hd[i].ps[j].Mencgas  += hd[i].ps[k].Mgas;
		hd[i].ps[j].Mencdark += hd[i].ps[k].Mdark;
		hd[i].ps[j].Mencstar += hd[i].ps[k].Mstar;
		}
	    for (k = 0; k < 3; k++) {
		hd[i].ps[j].veltot[k]  /= hd[i].ps[j].Ntot;
		hd[i].ps[j].velgas[k]  /= hd[i].ps[j].Ngas;
		hd[i].ps[j].veldark[k] /= hd[i].ps[j].Ndark;
		hd[i].ps[j].velstar[k] /= hd[i].ps[j].Nstar;
		}
	    for (k = 0; k < 6; k++) {
		hd[i].ps[j].vel2tot[k]  /= hd[i].ps[j].Ntot;
		hd[i].ps[j].vel2gas[k]  /= hd[i].ps[j].Ngas;
		hd[i].ps[j].vel2dark[k] /= hd[i].ps[j].Ndark;
		hd[i].ps[j].vel2star[k] /= hd[i].ps[j].Nstar;
		}
	    for (k = 0; k < 3; k++) {
		hd[i].ps[j].vel2tot[k]  -= hd[i].ps[j].veltot[k]*hd[i].ps[j].veltot[k];
		hd[i].ps[j].vel2gas[k]  -= hd[i].ps[j].velgas[k]*hd[i].ps[j].velgas[k];
		hd[i].ps[j].vel2dark[k] -= hd[i].ps[j].veldark[k]*hd[i].ps[j].veldark[k];
		hd[i].ps[j].vel2star[k] -= hd[i].ps[j].velstar[k]*hd[i].ps[j].velstar[k];
		}
	    hd[i].ps[j].vel2tot[3]  -= hd[i].ps[j].veltot[0]*hd[i].ps[j].veltot[1];
	    hd[i].ps[j].vel2tot[4]  -= hd[i].ps[j].veltot[0]*hd[i].ps[j].veltot[2];
	    hd[i].ps[j].vel2tot[5]  -= hd[i].ps[j].veltot[1]*hd[i].ps[j].veltot[2];
	    hd[i].ps[j].vel2gas[3]  -= hd[i].ps[j].velgas[0]*hd[i].ps[j].velgas[1];
	    hd[i].ps[j].vel2gas[4]  -= hd[i].ps[j].velgas[0]*hd[i].ps[j].velgas[2];
	    hd[i].ps[j].vel2gas[5]  -= hd[i].ps[j].velgas[1]*hd[i].ps[j].velgas[2];
	    hd[i].ps[j].vel2dark[3] -= hd[i].ps[j].veldark[0]*hd[i].ps[j].veldark[1];
	    hd[i].ps[j].vel2dark[4] -= hd[i].ps[j].veldark[0]*hd[i].ps[j].veldark[2];
	    hd[i].ps[j].vel2dark[5] -= hd[i].ps[j].veldark[1]*hd[i].ps[j].veldark[2];
	    hd[i].ps[j].vel2star[3] -= hd[i].ps[j].velstar[0]*hd[i].ps[j].velstar[1];
	    hd[i].ps[j].vel2star[4] -= hd[i].ps[j].velstar[0]*hd[i].ps[j].velstar[2];
	    hd[i].ps[j].vel2star[5] -= hd[i].ps[j].velstar[1]*hd[i].ps[j].velstar[2];
	    }
	
	/*
	** Calculate rbg, Mbg, rcrit, Mcrit by going from outside in
	*/
	
	for (j = gi.Nbin; j > 0; j--) {
	    radius[0] = hd[i].ps[j-1].ro;
	    radius[1] = hd[i].ps[j].ro;
	    rhoenc[0] = 3*hd[i].ps[j-1].Menctot/(4*M_PI*hd[i].ps[j-1].ro*hd[i].ps[j-1].ro*hd[i].ps[j-1].ro);
	    rhoenc[1] = 3*hd[i].ps[j].Menctot/(4*M_PI*hd[i].ps[j].ro*hd[i].ps[j].ro*hd[i].ps[j].ro);
	    Menc[0] = hd[i].ps[j-1].Menctot;
	    Menc[1] = hd[i].ps[j].Menctot;
	    if ((rhoenc[0] >= gi.rhoencbg) && (rhoenc[1] < gi.rhoencbg) && (hd[i].rbg == 0)) {
		if (gi.gridtype == 0) {
		    m = (radius[1]-radius[0])/(rhoenc[1]-rhoenc[0]);
		    d = gi.rhoencbg-rhoenc[0];
		    hd[i].rbg = radius[0] + m*d;
		    m = (Menc[1]-Menc[2])/(radius[1]-radius[0]);
		    d = hd[i].rbg-radius[0];
		    hd[i].Mbg = Menc[0] + m*d;
		    }
		else if (gi.gridtype == 1) {
		    m = (log(radius[1])-log(radius[0]))/(log(rhoenc[1])-log(rhoenc[0]));
		    d = log(gi.rhoencbg)-log(rhoenc[0]);
		    hd[i].rbg = exp(log(radius[0]) + m*d);
		    m = (log(Menc[1])-log(Menc[0]))/(log(radius[1])-log(radius[0]));
		    d = log(hd[i].rbg)-log(radius[0]);
		    hd[i].Mbg = exp(log(Menc[0]) + m*d);
		    }
		}
	    if ((rhoenc[0] >= gi.rhoenccrit) && (rhoenc[1] < gi.rhoenccrit) && (hd[i].rcrit == 0)) {
		if (gi.gridtype == 0) {
		    m = (radius[1]-radius[0])/(rhoenc[1]-rhoenc[0]);
		    d = gi.rhoenccrit-rhoenc[0];
		    hd[i].rcrit = radius[0] + m*d;
		    m = (Menc[1]-Menc[0])/(radius[1]-radius[0]);
		    d = hd[i].rcrit-radius[0];
		    hd[i].Mcrit = Menc[0] + m*d;
		    }
		else if (gi.gridtype == 1) {
		    m = (log(radius[1])-log(radius[0]))/(log(rhoenc[1])-log(rhoenc[0]));
		    d = log(gi.rhoenccrit)-log(rhoenc[0]);
		    hd[i].rcrit = exp(log(radius[0]) + m*d);
		    m = (log(Menc[1])-log(Menc[0]))/(log(radius[1])-log(radius[0]));
		    d = log(hd[i].rcrit)-log(radius[0]);
		    hd[i].Mcrit = exp(log(Menc[0]) + m*d);
		    }
		}
	    }

	/*
	** Calcualte vcmax, rvcmax, Mrvcmax by going from inside out
	*/
	
	for (j = 2; j < (gi.Nbin+1); j++) {
	    radius[0] = hd[i].ps[j-1].rm;
	    radius[1] = hd[i].ps[j].rm;
	    rhoenc[0] = 4*M_PI*(hd[i].ps[j-1].Mtot/hd[i].ps[j-1].vol)*hd[i].ps[j-1].rm*hd[i].ps[j-1].rm*hd[i].ps[j-1].rm;
	    rhoenc[1] = 4*M_PI*(hd[i].ps[j].Mtot/hd[i].ps[j].vol)*hd[i].ps[j].rm*hd[i].ps[j].rm*hd[i].ps[j].rm;
	    if (gi.gridtype == 0) {
		m = (hd[i].ps[j-1].Menctot-hd[i].ps[j-2].Menctot)/(hd[i].ps[j-1].ro-hd[i].ps[j-2].ro);
		d = radius[0]-hd[i].ps[j-2].ro;
		Menc[0] = hd[i].ps[j-2].Menctot + m*d;
		m = (hd[i].ps[j].Menctot-hd[i].ps[j-1].Menctot)/(hd[i].ps[j].ro-hd[i].ps[j-1].ro);
		d = radius[1]-hd[i].ps[j-1].ro;
		Menc[1] = hd[i].ps[j-1].Menctot + m*d;
		}
	    else if (gi.gridtype == 1) {
		m = (log(hd[i].ps[j-1].Menctot)-log(hd[i].ps[j-2].Menctot))/(log(hd[i].ps[j-1].ro)-log(hd[i].ps[j-2].ro));
		d = log(radius[0])-log(hd[i].ps[j-2].ro);
		Menc[0] = exp(log(hd[i].ps[j-2].Menctot) + m*d);
		m = (log(hd[i].ps[j].Menctot)-log(hd[i].ps[j-1].Menctot))/(log(hd[i].ps[j].ro)-log(hd[i].ps[j-1].ro));
		d = log(radius[1])-log(hd[i].ps[j-1].ro);
		Menc[1] = exp(log(hd[i].ps[j-1].Menctot) + m*d);
		}
	    if ((rhoenc[0] >= Menc[0]) && (rhoenc[1] < Menc[1]) && (hd[i].rvcmax == 0)) {
		rhoenc[0] -= Menc[0];
		rhoenc[1] -= Menc[1];
		if (gi.gridtype == 0) {
		    m = (radius[1]-radius[0])/(rhoenc[1]-rhoenc[0]);
		    d = 0 - rhoenc[0];
		    hd[i].rvcmax = radius[0] + m*d;
		    m = (Menc[1]-Menc[0])/(radius[1]-radius[0]);
		    d = hd[i].rvcmax-radius[0];
		    hd[i].Mrvcmax = Menc[0] + m*d;
		    hd[i].vcmax = sqrt(gi.G*hd[i].Mrvcmax/hd[i].rvcmax);
		    }
		else if (gi.gridtype == 1) {
		    m = (log(radius[1])-log(radius[0]))/(rhoenc[1]-rhoenc[0]);
		    d = 0 - rhoenc[0];
		    hd[i].rvcmax = exp(log(radius[0]) + m*d);
		    m = (log(Menc[1])-log(Menc[0]))/(log(radius[1])-log(radius[0]));
		    d = log(hd[i].rvcmax)-log(radius[0]);
		    hd[i].Mrvcmax = exp(log(Menc[0]) + m*d);
		    hd[i].vcmax = sqrt(gi.G*hd[i].Mrvcmax/hd[i].rvcmax);
		    }
		}
	    }
	}
    }

void write_output(HD *hd, GI gi) {

    int i, j, k;
    char statisticsfilename[256], profilesfilename[256];
    FILE *statisticsfile, *profilesfile;

    sprintf(statisticsfilename,"%s.statistics",gi.outputname);
    statisticsfile = fopen(statisticsfilename,"w");
    assert(statisticsfile != NULL);
    fprintf(statisticsfile,"# GID rbg Mbg rcrit Mcrit rvcmax Mrvcmax vcmax / total 8 columns\n");
    for (i = 0; i < gi.NHalo; i++) {
	fprintf(statisticsfile,"%d %.6e %.6e %.6e %.6e %.6e %.6e %.6e\n",hd[i].ID,hd[i].rbg,hd[i].Mbg,hd[i].rcrit,hd[i].Mcrit,hd[i].rvcmax,hd[i].Mrvcmax,hd[i].vcmax);
	}
    fclose(statisticsfile);

    sprintf(profilesfilename,"%s.profiles",gi.outputname);
    profilesfile = fopen(profilesfilename,"w");
    assert(profilesfile != NULL);
    fprintf(profilesfile,"# GID ri rm ro vol Mtot Mgas Mdark Mstar Menctot Mencgas Mencdark Mencstar rhoMtot rhoMgas rhoMdark rhoMstar Ntot Ngas Ndark Nstar Nenctot Nencgas Nencdark Nencstar rhoNtot rhoNgas rhoNdark rhoNstar veltotx veltoty veltotz vel2totxx vel2totyy vel2totzz vel2totxy vel2totxz vel2totyz velgasx velgasy velgasz vel2gasxx vel2gasyy vel2gaszz vel2gasxy vel2gasxz vel2gasyz veldarkx veldarky veldarkz vel2darkxx vel2darkyy vel2darkzz vel2darkxy vel2darkxz vel2darkyz velstarx velstary velstarz vel2starxx vel2staryy vel2starzz vel2starxy vel2starxz vel2staryz Ltotx Ltoty Ltotz Lgasx Lgasy Lgasz Ldarkx Ldarky Ldarkz Lstarx Lstary Lstarz / total 77 columns\n");
    
    for (i = 0; i < gi.NHalo; i++) {
	for (j = 0; j < (gi.Nbin+1); j++) {
	    fprintf(profilesfile,"%d ",hd[i].ID);
	    fprintf(profilesfile,"%.6e %.6e %.6e %.6e ",hd[i].ps[j].ri,hd[i].ps[j].rm,hd[i].ps[j].ro,hd[i].ps[j].vol);
	    fprintf(profilesfile,"%.6e %.6e %.6e %.6e ",hd[i].ps[j].Mtot,hd[i].ps[j].Mgas,hd[i].ps[j].Mdark,hd[i].ps[j].Mstar);
	    fprintf(profilesfile,"%.6e %.6e %.6e %.6e ",hd[i].ps[j].Menctot,hd[i].ps[j].Mencgas,hd[i].ps[j].Mencdark,hd[i].ps[j].Mencstar);
	    fprintf(profilesfile,"%.6e %.6e %.6e %.6e ",hd[i].ps[j].Mtot/hd[i].ps[j].vol,hd[i].ps[j].Mgas/hd[i].ps[j].vol,hd[i].ps[j].Mdark/hd[i].ps[j].vol,hd[i].ps[j].Mstar/hd[i].ps[j].vol);
	    fprintf(profilesfile,"%d %d %d %d ",hd[i].ps[j].Ntot,hd[i].ps[j].Ngas,hd[i].ps[j].Ndark,hd[i].ps[j].Nstar);
	    fprintf(profilesfile,"%d %d %d %d ",hd[i].ps[j].Nenctot,hd[i].ps[j].Nencgas,hd[i].ps[j].Nencdark,hd[i].ps[j].Nencstar);
	    fprintf(profilesfile,"%.6e %.6e %.6e %.6e ",hd[i].ps[j].Ntot/hd[i].ps[j].vol,hd[i].ps[j].Ngas/hd[i].ps[j].vol,hd[i].ps[j].Ndark/hd[i].ps[j].vol,hd[i].ps[j].Nstar/hd[i].ps[j].vol);
	    for (k = 0; k < 3; k++) fprintf(profilesfile,"%.6e ",hd[i].ps[j].veltot[k]);
	    for (k = 0; k < 6; k++) fprintf(profilesfile,"%.6e ",hd[i].ps[j].vel2tot[k]);
	    for (k = 0; k < 3; k++) fprintf(profilesfile,"%.6e ",hd[i].ps[j].velgas[k]);
	    for (k = 0; k < 6; k++) fprintf(profilesfile,"%.6e ",hd[i].ps[j].vel2gas[k]);
	    for (k = 0; k < 3; k++) fprintf(profilesfile,"%.6e ",hd[i].ps[j].veldark[k]);
	    for (k = 0; k < 6; k++) fprintf(profilesfile,"%.6e ",hd[i].ps[j].vel2dark[k]);
	    for (k = 0; k < 3; k++) fprintf(profilesfile,"%.6e ",hd[i].ps[j].velstar[k]);
	    for (k = 0; k < 6; k++) fprintf(profilesfile,"%.6e ",hd[i].ps[j].vel2star[k]);
	    for (k = 0; k < 3; k++) fprintf(profilesfile,"%.6e ",hd[i].ps[j].Ltot[k]);
	    for (k = 0; k < 3; k++) fprintf(profilesfile,"%.6e ",hd[i].ps[j].Lgas[k]);
	    for (k = 0; k < 3; k++) fprintf(profilesfile,"%.6e ",hd[i].ps[j].Ldark[k]);
	    for (k = 0; k < 3; k++) fprintf(profilesfile,"%.6e ",hd[i].ps[j].Lstar[k]);
	    fprintf(profilesfile,"\n");
	    }
	}
    fclose(profilesfile);
    }
