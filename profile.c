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

#define ART_LEVEL_ARRAY_LENGTH 10
#define ART_BANNER_LENGTH 45

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
    int Nbin;
    double rcentre[3];
    double vcentre[3];
    double rbg, Mbg;
    double rcrit, Mcrit;
    double rvcmax, Mrvcmax, vcmax;
    double rmin, rmax;
    PS *ps;
    } HALO_DATA;

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
    int Nbin;
    int NHalo, Ngas, Ndark, Nstar;
    int NParticlesPerBlock, NParticlesInBlock, NBlock;
    double rhoencbg, rhoenccrit;
    double rmin, rmax;
    double E, Deltabg, Deltacrit;
    double bc[6];
    double binfactor;
    char HaloCatalogueFileName[256];
    char OutputName[256];
    } GI;


// => goes to iof
typedef struct cosmological_parameters {

    double ascale0;
    double OmegaM0;
    double OmegaDM0; 
    double OmegaB0;
    double OmegaL0;
    double OmegaK0;
    double OmegaR0;
    double Hubble0;
    double ascale;
    double OmegaM;
    double OmegaDM; 
    double OmegaB;
    double OmegaL;
    double OmegaK;
    double OmegaR;
    double Hubble;
    } COSMOLOGICAL_PARAMETERS;

// => goes to iof
typedef struct unit_system {

    double GNewton;
    double rhocrit0;
    double rhocrit;
    double LBox;
    } UNIT_SYSTEM;


// => goes to iof
typedef struct art_data {

    int doswap;
    int massfromdata;
    int Lmaxdark, Nrec;
    int headerincludesstars;
    int Ngas, Ndark, Nstar;
    int Nreadgas, Nreaddark, Nreadstar;
    int Ndarklevel[ART_LEVEL_ARRAY_LENGTH];
    double shift;
    double toplevelmassdark, toplevelsoftdark, refinementstepdark;
    double massdark[ART_LEVEL_ARRAY_LENGTH];
    double softdark[ART_LEVEL_ARRAY_LENGTH];
    char HeaderFileName[256], GasFileName[256], DarkFileName[256], StarFileName[256];
    char Banner[ART_BANNER_LENGTH];
    ART_HEADER ah;
    FILE *HeaderFile;
    FILE *PosXFile, *PosYFile, *PosZFile;
    FILE *VelXFile, *VelYFile, *VelZFile;
    } ART_DATA;

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


// => goes to iof
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

void read_art_header(ART_DATA *ad) {

    int i, L;
    int header, trailer;

    ad->doswap = 0;
    ad->HeaderFile = fopen(ad->HeaderFileName,"r");
    assert(ad->HeaderFile != NULL);
    assert(fread(&header,sizeof(int),1,ad->HeaderFile) == 1);
    if (header != ART_BANNER_LENGTH+sizeof(ART_HEADER)) {
	ad->doswap = 1;
	flip_4byte(&header,sizeof(int),1);
	}
    assert(header == ART_BANNER_LENGTH+sizeof(ART_HEADER));
    assert(fread(ad->Banner,sizeof(char),ART_BANNER_LENGTH,ad->HeaderFile) == ART_BANNER_LENGTH);
    assert(fread(&ad->ah,sizeof(ART_HEADER),1,ad->HeaderFile) == 1);
    if (ad->doswap) flip_4byte(&ad->ah,sizeof(ART_HEADER),1);
    assert(fread(&trailer,sizeof(int),1,ad->HeaderFile) == 1);
    if (ad->doswap) flip_4byte(&trailer,sizeof(int),1);
    assert(header == trailer);
    fclose(ad->HeaderFile);
    /*
    ** Set some derived quantities
    */
    ad->Ngas = 0;
    if (ad->headerincludesstars == 1) {
	ad->Lmaxdark = ad->ah.Nspecies-2;
	ad->Ndark = ad->ah.num[ad->Lmaxdark];
	ad->Nstar = ad->ah.num[ad->ah.Nspecies-1];
	}
    else {
	ad->Lmaxdark = ad->ah.Nspecies-1;
	ad->Ndark = ad->ah.num[ad->Lmaxdark];
	ad->Nstar = 0;
	}
    ad->Nrec = ad->ah.Nrow*ad->ah.Nrow;
    if (ad->toplevelmassdark == -1) {
	ad->massfromdata = 1;
	ad->toplevelmassdark = ad->ah.mass[ad->Lmaxdark];
	}
    assert(ad->toplevelsoftdark >= 0);
    assert(ad->toplevelmassdark >= 0);
    for (i = 0; i <= ad->Lmaxdark; i++) {
	L = ad->Lmaxdark-i;
	if (i == 0) {
	    ad->Ndarklevel[L] = ad->ah.num[i];
	    }
	else {
	    ad->Ndarklevel[L] = ad->ah.num[i]-ad->ah.num[i-1];
	    }
	if (ad->massfromdata == 1) {
	    ad->massdark[L] = ad->ah.mass[i];
	    }
	}
    }

void usage(void);
void read_halocatalogue_6DFOF(GI (*), HALO_DATA (**));
void read_art_record_dark(ART_DATA, DARK_PARTICLE (*));
void put_dp_in_bins(HALO_DATA (*), DARK_PARTICLE (*), UNIT_SYSTEM, GI);
void calculate_halo_properties(HALO_DATA (*), COSMOLOGICAL_PARAMETERS, UNIT_SYSTEM, GI);
void write_output(HALO_DATA (*), GI);

int main(int argc, char **argv) {

    int i, j;
    GI gi;
    COSMOLOGICAL_PARAMETERS cp;
    UNIT_SYSTEM us;
    TIPSY_HEADER th;
    GAS_PARTICLE *gp = NULL;
    DARK_PARTICLE *dp = NULL;
    STAR_PARTICLE *sp = NULL;
    GAS_PARTICLE_DPP *gpdpp = NULL;
    DARK_PARTICLE_DPP *dpdpp = NULL;
    STAR_PARTICLE_DPP *spdpp = NULL;
    HALO_DATA *hd = NULL;
    ART_DATA ad;
    FILE *halocataloguefile;
    XDR xdrs;

    /*
    ** Set some default values
    */

    cp.ascale0 = 1;
    cp.OmegaM0 = 0.3;
    cp.OmegaDM0 = 0.26;
    cp.OmegaB0 = 0.04;
    cp.OmegaL0 = 0.7;
    cp.OmegaK0 = 0;
    cp.OmegaR0 = 0;
    cp.Hubble0 = sqrt(8.0*M_PI/3.0);

    us.GNewton = 1;
    us.rhocrit0 = 1;
    us.LBox = 1;

    gi.Deltabg = 200;
    gi.Deltacrit = 0;
    gi.rmin = 0;
    gi.rmax = 1;
    gi.Nbin = 0;
    gi.binfactor = 5;
    gi.NHalo = 0;

    ad.doswap = 0;
    ad.massfromdata = 0;
    ad.headerincludesstars = 0;
    ad.Ngas = 0;
    ad.Ndark = 0;
    ad.Nstar = 0;
    ad.Nreadgas = 0;
    ad.Nreaddark = 0;
    ad.Nreadstar = 0;
    ad.refinementstepdark = 2;
    ad.toplevelmassdark = -1;
    ad.toplevelsoftdark = 0;
    ad.shift = 0;
    for (i = 0; i < ART_LEVEL_ARRAY_LENGTH; i++) {
	ad.Ndarklevel[i] = 0;
	ad.massdark[i] = 0;
	ad.softdark[i] = 0;
	}

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
            cp.OmegaM0 = atof(argv[i]);
            i++;
            }
        else if (strcmp(argv[i],"-OmegaK0") == 0) {
            i++;
            if (i >= argc) usage();
            cp.OmegaK0 = atof(argv[i]);
            i++;
            }
        else if (strcmp(argv[i],"-OmegaL0") == 0) {
            i++;
            if (i >= argc) usage();
            cp.OmegaL0 = atof(argv[i]);
            i++;
            }
        else if (strcmp(argv[i],"-LBox") == 0) {
            i++;
            if (i >= argc) usage();
            us.LBox = atof(argv[i]);
            i++;
            }
        else if (strcmp(argv[i],"-rhocrit0") == 0) {
            i++;
            if (i >= argc) usage();
            us.rhocrit0 = atof(argv[i]);
            i++;
            }
        else if (strcmp(argv[i],"-G") == 0) {
            i++;
            if (i >= argc) usage();
            us.GNewton = atof(argv[i]);
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
            gi.rmaxfromhalocatalogue = 0;
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
            strcpy(gi.HaloCatalogueFileName,argv[i]);
            i++;
            }
        else if (strcmp(argv[i],"-output") == 0) {
            i++;
            if (i >= argc) usage();
            strcpy(gi.OutputName,argv[i]);
            i++;
            }
        else if (strcmp(argv[i],"-headerfile") == 0) {
            i++;
            if (i >= argc) usage();
            strcpy(ad.HeaderFileName,argv[i]);
            i++;
            }
        else if (strcmp(argv[i],"-gasfile") == 0) {
            i++;
            if (i >= argc) usage();
            strcpy(ad.GasFileName,argv[i]);
            i++;
            }
        else if (strcmp(argv[i],"-darkfile") == 0) {
            i++;
            if (i >= argc) usage();
            strcpy(ad.DarkFileName,argv[i]);
            i++;
            }
        else if (strcmp(argv[i],"-starfile") == 0) {
            i++;
            if (i >= argc) usage();
            strcpy(ad.StarFileName,argv[i]);
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
	cp.ascale = th.time;
	gi.Ngas = th.ngas;
	gi.Ndark = th.ndark;
	gi.Nstar = th.nstar;
	gi.NParticlesPerBlock = 1000000;
	gi.bc[0] = -0.5*us.LBox;
	gi.bc[1] = -0.5*us.LBox;
	gi.bc[2] = -0.5*us.LBox;
	gi.bc[3] = 0.5*us.LBox;
	gi.bc[4] = 0.5*us.LBox;
	gi.bc[5] = 0.5*us.LBox;
	}
    else if (gi.dataformat == 1) {
	read_art_header(&ad);
	cp.ascale = ad.ah.aunin;
	gi.Ngas = ad.Ngas;
	gi.Ndark = ad.Ndark;
	gi.Nstar = ad.Nstar;
	cp.OmegaM0 = ad.ah.OmM0;
	cp.OmegaB0 = ad.ah.OmB0;
	cp.OmegaDM0 = cp.OmegaM0 - cp.OmegaB0;
	cp.OmegaL0 = ad.ah.OmL0;
	cp.OmegaK0 = ad.ah.OmK0;
	gi.NParticlesPerBlock = ad.Nrec;
	gi.bc[0] = 0;
	gi.bc[1] = 0;
	gi.bc[2] = 0;
	gi.bc[3] = us.LBox;
	gi.bc[4] = us.LBox;
	gi.bc[5] = us.LBox;
	}
    else {
	fprintf(stderr,"Not supported format!\n");
	exit(1);
	}

    fprintf(stderr,"a = %g doswap %d Nrec %d %d B %s Ndark %d\n",cp.ascale,ad.doswap,ad.Nrec,ad.ah.Nrow,ad.Banner,ad.Ndark);

    for (i = 0; i < 10; i++) {
	fprintf(stderr,"i %d Ndarklvel %d massdark %g softdark %g num %d\n",i,ad.Ndarklevel[i],ad.massdark[i],ad.softdark[i],ad.ah.num[i]);
	}

//    exit(1);

    /*
    ** Calculate cosmology relevant stuff => write function for this
    ** Densities in comoving coordinates 
    */

    gi.rhoencbg = gi.Deltabg*cp.OmegaM0*us.rhocrit0;
    gi.E = Ecosmo(cp.ascale,cp.OmegaM0,cp.OmegaL0,cp.OmegaK0);
    cp.OmegaM = cp.OmegaM0/(pow(cp.ascale,3)*gi.E*gi.E);
    us.rhocrit = us.rhocrit0*gi.E*gi.E;
    if (gi.Deltacrit == 0) {
	if (cp.OmegaK0 == 0) {
	    gi.Deltacrit = 178*pow(cp.OmegaM,0.45);
	    }
	else if (cp.OmegaL0 == 0) {
	    gi.Deltacrit = 178*pow(cp.OmegaM,0.3);
	    }
	}
    gi.rhoenccrit = gi.Deltacrit*us.rhocrit*pow(cp.ascale,3);

    fprintf(stderr,"a = %g doswap %d \n",cp.ascale,ad.doswap);

    /*
    ** Read halo catalogue
    */

    if (gi.halocatalogueformat == 0) {
//	read_halocatalogue_generic(gi,&hd);
	}
    else if (gi.halocatalogueformat == 1) {
	read_halocatalogue_6DFOF(&gi,&hd);
	}


    fprintf(stderr,"After halocatalogue\n");
    fprintf(stderr,"file> %s NHalo %d\n",gi.HaloCatalogueFileName,gi.NHalo);
    for (i = 0; i < gi.NHalo; i++) {
	fprintf(stderr,"%d %g %g %g %g %g %g %g %g %d\n",hd[i].ID,hd[i].rcentre[0],hd[i].rcentre[1],hd[i].rcentre[2],hd[i].vcentre[0],hd[i].vcentre[1],hd[i].vcentre[2],hd[i].rmin,hd[i].rmax,hd[i].Nbin);
	}

    exit(1);

    /*
    ** Harvest data
    */
    if (gi.dataformat == 0) {
	assert(gi.Ngas == 0);
	assert(gi.Nstar == 0);
	assert(gi.positionprecision == 0);
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
		j = gi.NParticlesInBlock-1;
		fprintf(stderr,"%d %d %g %g %g\n",i,j,dp[j].pos[0],dp[j].pos[1],dp[j].pos[2]);

		put_dp_in_bins(hd,dp,us,gi);
		}
	    fprintf(stderr,"Ndark %d = %d\n",gi.Ndark,(gi.NBlock-1)*gi.NParticlesPerBlock+gi.NParticlesInBlock);
	    free(dp);
	    }
	}
    else if (gi.dataformat == 1) {
	/*
	** Dark Matter
	*/
	ad.PosXFile = fopen(ad.DarkFileName,"r");
	assert(ad.PosXFile != NULL);
	ad.PosYFile = fopen(ad.DarkFileName,"r");
	assert(ad.PosYFile != NULL);
	ad.PosZFile = fopen(ad.DarkFileName,"r");
	assert(ad.PosZFile != NULL);
	ad.VelXFile = fopen(ad.DarkFileName,"r");
	assert(ad.VelXFile != NULL);
	ad.VelYFile = fopen(ad.DarkFileName,"r");
	assert(ad.VelYFile != NULL);
	ad.VelZFile = fopen(ad.DarkFileName,"r");
	assert(ad.VelZFile != NULL);
	dp = realloc(dp,gi.NParticlesPerBlock*sizeof(DARK_PARTICLE));
	gi.NBlock = (gi.Ndark+gi.NParticlesPerBlock-1)/gi.NParticlesPerBlock;
	for (i = 0; i < gi.NBlock; i++) {
	    if (i == gi.NBlock-1) {
		gi.NParticlesInBlock = gi.Ndark-((gi.NBlock-1)*gi.NParticlesPerBlock);
		}
	    else {
		gi.NParticlesInBlock = gi.NParticlesPerBlock;
		}
	    read_art_record_dark(ad,dp);
	    j = gi.NParticlesInBlock-1;
	    fprintf(stderr,"%d %d %g %g %g\n",i,j,dp[j].pos[0],dp[j].pos[1],dp[j].pos[2]);
	    put_dp_in_bins(hd,dp,us,gi);
	    }
	fprintf(stderr,"Ndark %d = %d\n",gi.Ndark,(gi.NBlock-1)*gi.NParticlesPerBlock+gi.NParticlesInBlock);
	fclose(ad.PosXFile);
	fclose(ad.PosYFile);
	fclose(ad.PosZFile);
	fclose(ad.VelXFile);
	fclose(ad.VelYFile);
	fclose(ad.VelZFile);
	/*
	** Stars
	*/

	}

    /*
    ** Calculate halo properties
    */
 
    calculate_halo_properties(hd,cp,us,gi); 

    /*
    ** Write output
    */

    write_output(hd,gi);

    /*
    ** Some more output if desired
    */


    if (gi.verboselevel >= 1) {
        fprintf(stderr,"Used values:\n");
        fprintf(stderr,"Delta_bg    : %.6e\n",gi.Deltabg);
        fprintf(stderr,"Delta_crit  : %.6e\n",gi.Deltacrit);
	fprintf(stderr,"rhoenc_bg   : %.6e MU LU^{-3} (comoving)\n",gi.rhoencbg);
	fprintf(stderr,"rhoenc_crit : %.6e MU LU^{-3} (comoving)\n",gi.rhoenccrit);
        fprintf(stderr,"rhocrit0    : %.6e MU LU^{-3}\n",us.rhocrit0);
        fprintf(stderr,"OmegaM0     : %.6e\n",cp.OmegaM0);
        fprintf(stderr,"OmegaL0     : %.6e\n",cp.OmegaL0);
        fprintf(stderr,"OmegaK0     : %.6e\n",cp.OmegaK0);
        fprintf(stderr,"LBox        : %.6e\n",us.LBox);
        fprintf(stderr,"G           : %.6e LU^3 MU^{-1} TU^{-1}\n",us.GNewton);
	fprintf(stderr,"a           : %.6e\n",cp.ascale);
	fprintf(stderr,"E           : %.6e\n",gi.E);
        fprintf(stderr,"binfactor   : %.6e\n",gi.binfactor);
        }

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

void read_halocatalogue_6DFOF(GI *gi, HALO_DATA **hdin) {

    int SizeHaloData = 10000;
    int i, j, k, ID, N, idummy, NHaloRead;
    float fdummy;
    double DarkMass, radius1, radius2, vd1D;
    double rxcom, rycom, rzcom, rxpotmin, rypotmin, rzpotmin, rxdenmax, rydenmax, rzdenmax, vx, vy, vz;
    double dr;
    HALO_DATA *hd;
    FILE *HaloCatalogueFile = NULL;

    fprintf(stderr,"pointer1 %p\n",HaloCatalogueFile);

    HaloCatalogueFile = fopen(gi->HaloCatalogueFileName,"r");
    assert(HaloCatalogueFile != NULL);

    fprintf(stderr,"pointer1 %p\n",HaloCatalogueFile);

    hd = *hdin;
    assert(gi->Nbin > 0);
    assert(gi->rmin >= 0);
    assert(gi->rmax >= 0);
    assert(gi->rmax > gi->rmin);
    dr = 0;
    NHaloRead = 0;
    hd = realloc(hd,SizeHaloData*sizeof(HALO_DATA));
    assert(hd != NULL);
    for (i = 0; i < SizeHaloData; i++) {
	hd[i].Nbin = gi->Nbin;
	hd[i].ps = realloc(hd[i].ps,(hd[i].Nbin+1)*sizeof(PS));
	assert(hd[i].ps != NULL);
	}
    while (NHaloRead < 2) {
	fscanf(HaloCatalogueFile,"%i",&idummy); ID = idummy;
	fscanf(HaloCatalogueFile,"%i",&idummy); N = idummy;
	fscanf(HaloCatalogueFile,"%g",&fdummy); DarkMass = fdummy;
	fscanf(HaloCatalogueFile,"%g",&fdummy);
	fscanf(HaloCatalogueFile,"%g",&fdummy);
	fscanf(HaloCatalogueFile,"%g",&fdummy);
	fscanf(HaloCatalogueFile,"%g",&fdummy); radius1 = fdummy; /* (sum_j rmax[j]-rmin[j])/6 */
	fscanf(HaloCatalogueFile,"%g",&fdummy); radius2 = fdummy; /* dispersion in coordinates */
	fscanf(HaloCatalogueFile,"%g",&fdummy); vd1D = fdummy;
	fscanf(HaloCatalogueFile,"%g",&fdummy);
	fscanf(HaloCatalogueFile,"%g",&fdummy);
	fscanf(HaloCatalogueFile,"%g",&fdummy);
	fscanf(HaloCatalogueFile,"%g",&fdummy); rxcom = put_in_box(fdummy,gi->bc[0],gi->bc[3]);
	fscanf(HaloCatalogueFile,"%g",&fdummy); rycom = put_in_box(fdummy,gi->bc[1],gi->bc[4]);
	fscanf(HaloCatalogueFile,"%g",&fdummy); rzcom = put_in_box(fdummy,gi->bc[2],gi->bc[5]);
	fscanf(HaloCatalogueFile,"%g",&fdummy); rxpotmin = put_in_box(fdummy,gi->bc[0],gi->bc[3]);
	fscanf(HaloCatalogueFile,"%g",&fdummy); rypotmin = put_in_box(fdummy,gi->bc[1],gi->bc[4]);
	fscanf(HaloCatalogueFile,"%g",&fdummy); rzpotmin = put_in_box(fdummy,gi->bc[2],gi->bc[5]);
	fscanf(HaloCatalogueFile,"%g",&fdummy); rxdenmax = put_in_box(fdummy,gi->bc[0],gi->bc[3]);
	fscanf(HaloCatalogueFile,"%g",&fdummy); rydenmax = put_in_box(fdummy,gi->bc[1],gi->bc[4]);
	fscanf(HaloCatalogueFile,"%g",&fdummy); rzdenmax = put_in_box(fdummy,gi->bc[2],gi->bc[5]);
	fscanf(HaloCatalogueFile,"%g",&fdummy); vx = fdummy;
	fscanf(HaloCatalogueFile,"%g",&fdummy); vy = fdummy;
	fscanf(HaloCatalogueFile,"%g",&fdummy); vz = fdummy;
	fscanf(HaloCatalogueFile,"%g",&fdummy);
	fscanf(HaloCatalogueFile,"%g",&fdummy);
	fscanf(HaloCatalogueFile,"%g",&fdummy);
	fscanf(HaloCatalogueFile,"%g",&fdummy);
	fscanf(HaloCatalogueFile,"%g",&fdummy);
	if (feof(HaloCatalogueFile)) break;
	NHaloRead++;
	if (SizeHaloData < NHaloRead){
	    SizeHaloData += 10000; 
	    hd = realloc(hd,SizeHaloData*sizeof(HALO_DATA));
	    assert(hd != NULL);
	    for (i = 0; i < SizeHaloData; i++) {
		hd[i].Nbin = gi->Nbin;
		hd[i].ps = realloc(hd[i].ps,(hd[i].Nbin+1)*sizeof(PS));
		assert(hd[i].ps != NULL);
		}
	    }
	i = NHaloRead-1;
	hd[i].ID = ID;
	if (gi->centretype == 0) {
	    hd[i].rcentre[0] = rxcom;
	    hd[i].rcentre[1] = rycom;
	    hd[i].rcentre[2] = rzcom;
	    }
	else if (gi->centretype == 1) {
	    hd[i].rcentre[0] = rxpotmin;
	    hd[i].rcentre[1] = rypotmin;
	    hd[i].rcentre[2] = rzpotmin;
	    }
	else if (gi->centretype == 2) {
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
	hd[i].rmin = gi->rmin;
	if (gi->rmaxfromhalocatalogue == 1) {
	    /*
	    ** Estimate maximum radius; assume isothermal sphere scaling 
	    */
	    hd[i].rmax = sqrt((3*DarkMass/(4*M_PI*radius2*radius2*radius2))/gi->rhoencbg)*radius2*gi->binfactor;
	    assert(hd[i].rmax > 0);
	    }
	else {
	    hd[i].rmax = gi->rmax;
	    }
	if (gi->gridtype == 0) {
	    dr = (hd[i].rmax-hd[i].rmin)/hd[i].Nbin;
	    }
	else if (gi->gridtype == 1) {
	    dr = (log(hd[i].rmax)-log(hd[i].rmin))/hd[i].Nbin;
	    }
	hd[i].ps[0].ri = 0;
	hd[i].ps[0].ro = hd[i].rmin;
	for (j = 1; j < (hd[i].Nbin+1); j++) {
	    if (gi->gridtype == 0) {
		hd[i].ps[j].ri = hd[i].rmin + (j-1)*dr;
		hd[i].ps[j].ro = hd[i].rmin + j*dr;
		}
	    else if (gi->gridtype == 1) {
		hd[i].ps[j].ri = exp(log(hd[i].rmin) + (j-1)*dr);
		hd[i].ps[j].ro = exp(log(hd[i].rmin) + j*dr);
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
    if (HaloCatalogueFile == NULL) fprintf(stderr,"Pointer is NULL\n");
    fprintf(stderr,"pointer2 %p\n",HaloCatalogueFile);
//    rewind(HaloCatalogueFile);
    fprintf(stderr,"pointer3 %p\n",HaloCatalogueFile);
    fclose(HaloCatalogueFile);
    *hdin = hd;
    gi->NHalo = NHaloRead;

    fprintf(stderr,"file> %s NHalo %d\n",gi->HaloCatalogueFileName,gi->NHalo);
    for (i = 0; i < gi->NHalo; i++) {
	fprintf(stderr,"%d %g %g %g %g %g %g %g %g %d\n",hd[i].ID,hd[i].rcentre[0],hd[i].rcentre[1],hd[i].rcentre[2],hd[i].vcentre[0],hd[i].vcentre[1],hd[i].vcentre[2],hd[i].rmin,hd[i].rmax,hd[i].Nbin);
	}

    }

void read_art_record_dark(ART_DATA ad, DARK_PARTICLE *dp) {

    int j, k;
    int L = 0;
    double rx, ry, rz;
    double vx, vy, vz;
    
    /*
    ** Get file pointers ready
    */
    
    assert(fseek(ad.PosYFile,1*ad.Nrec*sizeof(float),SEEK_CUR) == 0);
    assert(fseek(ad.PosZFile,2*ad.Nrec*sizeof(float),SEEK_CUR) == 0);
    assert(fseek(ad.VelXFile,3*ad.Nrec*sizeof(float),SEEK_CUR) == 0);
    assert(fseek(ad.VelYFile,4*ad.Nrec*sizeof(float),SEEK_CUR) == 0);
    assert(fseek(ad.VelZFile,5*ad.Nrec*sizeof(float),SEEK_CUR) == 0);

    /*
    ** Go through record
    */
    
    for (j = 0; j < ad.Nrec; j++) {

	assert(fread(&rx,sizeof(float),1,ad.PosXFile) == 1);
	assert(fread(&ry,sizeof(float),1,ad.PosYFile) == 1);
	assert(fread(&rz,sizeof(float),1,ad.PosZFile) == 1);
	assert(fread(&vx,sizeof(float),1,ad.VelXFile) == 1);
	assert(fread(&vy,sizeof(float),1,ad.VelYFile) == 1);
	assert(fread(&vz,sizeof(float),1,ad.VelZFile) == 1);

	if (ad.doswap) flip_4byte(&rx,sizeof(float),1);
	if (ad.doswap) flip_4byte(&ry,sizeof(float),1);
	if (ad.doswap) flip_4byte(&rz,sizeof(float),1);
	if (ad.doswap) flip_4byte(&vx,sizeof(float),1);
	if (ad.doswap) flip_4byte(&vy,sizeof(float),1);
	if (ad.doswap) flip_4byte(&vz,sizeof(float),1);

	/*
	** Determine current level
	*/

	for (k = ad.Lmaxdark; k >=0; k--) {
	    if (ad.ah.num[k] > ad.Nreaddark) L = ad.Lmaxdark-k;
	    }

	/*
	** Set particle properties
	*/

	dp[j].pos[0] = rx+ad.shift;
	dp[j].pos[1] = ry+ad.shift;
	dp[j].pos[2] = rz+ad.shift;
	dp[j].vel[0] = vx;
	dp[j].vel[1] = vy;
	dp[j].vel[2] = vz;
	dp[j].mass = ad.massdark[L];
	dp[j].eps = ad.softdark[L];
	dp[j].phi = 0;
	ad.Nreaddark++;
	}

    /*
    ** Move all pointers to the end of the record
    */

    assert(fseek(ad.PosXFile,5*ad.Nrec*sizeof(float),SEEK_CUR) == 0);
    assert(fseek(ad.PosYFile,4*ad.Nrec*sizeof(float),SEEK_CUR) == 0);
    assert(fseek(ad.PosZFile,3*ad.Nrec*sizeof(float),SEEK_CUR) == 0);
    assert(fseek(ad.VelXFile,2*ad.Nrec*sizeof(float),SEEK_CUR) == 0);
    assert(fseek(ad.VelYFile,1*ad.Nrec*sizeof(float),SEEK_CUR) == 0);
    }

void put_dp_in_bins(HALO_DATA *hd, DARK_PARTICLE *dp, UNIT_SYSTEM us, GI gi) {

    int i, j, k, l;
    double pos[3], vel[3], velproj[3];
    double erad[3], ephi[3], etheta[3];
    double r;

    for (i = 0; i < gi.NParticlesPerBlock; i++) {
	for (j = 0; j < gi.NHalo; j++) {
	    for (k = 0; k < 3; k++) {
		pos[k] = correct_position(hd[j].rcentre[k],dp[i].pos[k],us.LBox);
		pos[k] = pos[k]-hd[j].rcentre[k];
		}
	    r = sqrt(pos[0]*pos[0]+pos[1]*pos[1]+pos[2]*pos[2]);
	    if (r <= hd[j].ps[gi.Nbin].ro) {
		for (k = 0; k < 3; k++) {
		    vel[k] = dp[i].vel[k]-hd[j].vcentre[k];
		    }
		/*
		** Go through bins from outside inwards => larger bin volume further out
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

void calculate_halo_properties(HALO_DATA *hd, COSMOLOGICAL_PARAMETERS cp, UNIT_SYSTEM us, GI gi) {

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
		    hd[i].vcmax = sqrt(us.GNewton*hd[i].Mrvcmax/hd[i].rvcmax);
		    }
		else if (gi.gridtype == 1) {
		    m = (log(radius[1])-log(radius[0]))/(rhoenc[1]-rhoenc[0]);
		    d = 0 - rhoenc[0];
		    hd[i].rvcmax = exp(log(radius[0]) + m*d);
		    m = (log(Menc[1])-log(Menc[0]))/(log(radius[1])-log(radius[0]));
		    d = log(hd[i].rvcmax)-log(radius[0]);
		    hd[i].Mrvcmax = exp(log(Menc[0]) + m*d);
		    hd[i].vcmax = sqrt(us.GNewton*hd[i].Mrvcmax/hd[i].rvcmax);
		    }
		}
	    }
	}
    }

void write_output(HALO_DATA *hd, GI gi) {

    int i, j, k;
    char statisticsfilename[256], profilesfilename[256];
    FILE *statisticsfile, *profilesfile;

    sprintf(statisticsfilename,"%s.statistics",gi.OutputName);
    statisticsfile = fopen(statisticsfilename,"w");
    assert(statisticsfile != NULL);
    fprintf(statisticsfile,"#GID/1 rx/2 ry/3 rz/4 vx/5 vy/6 vz/7 rbg/8 Mbg/9 rcrit/10 Mcrit/11 rvcmax/12 Mrvcmax/13 vcmax/14\n");
    for (i = 0; i < gi.NHalo; i++) {
	fprintf(statisticsfile,"%d %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e\n",
		hd[i].ID,hd[i].rcentre[0],hd[i].rcentre[1],hd[i].rcentre[2],hd[i].vcentre[0],hd[i].vcentre[1],hd[i].vcentre[2],
		hd[i].rbg,hd[i].Mbg,hd[i].rcrit,hd[i].Mcrit,hd[i].rvcmax,hd[i].Mrvcmax,hd[i].vcmax);
	}
    fclose(statisticsfile);

    sprintf(profilesfilename,"%s.profiles",gi.OutputName);
    profilesfile = fopen(profilesfilename,"w");
    assert(profilesfile != NULL);
    fprintf(profilesfile,"#GID/1 ri/2 rm/3 ro/4 vol/5 Mtot/6 Mgas/7 Mdark/8 Mstar/9 Menctot/10 Mencgas/11 Mencdark/12 Mencstar/13 rhoMtot/14 rhoMgas/15 rhoMdark/16 rhoMstar/17 Ntot/18 Ngas/19 Ndark/20 Nstar/21 Nenctot/22 Nencgas/23 Nencdark/24 Nencstar/25 rhoNtot/26 rhoNgas/27 rhoNdark/28 rhoNstar/29 veltotx/30 veltoty/31 veltotz/32 vel2totxx/33 vel2totyy/34 vel2totzz/35 vel2totxy/36 vel2totxz/37 vel2totyz/38 velgasx/39 velgasy/40 velgasz/41 vel2gasxx/42 vel2gasyy/43 vel2gaszz/44 vel2gasxy/45 vel2gasxz/46 vel2gasyz/47 veldarkx/48 veldarky/49 veldarkz/50 vel2darkxx/51 vel2darkyy/52 vel2darkzz/53 vel2darkxy/54 vel2darkxz/55 vel2darkyz/56 velstarx/57 velstary/58 velstarz/59 vel2starxx/60 vel2staryy/61 vel2starzz/62 vel2starxy/63 vel2starxz/64 vel2staryz/65 Ltotx/66 Ltoty/67 Ltotz/68 Lgasx/69 Lgasy/70 Lgasz/71 Ldarkx/72 Ldarky/73 Ldarkz/74 Lstarx/75 Lstary/76 Lstarz/77\n");
    
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
