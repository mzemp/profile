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

typedef struct profilearray {

    int Ntot;
    int Ngas;
    int Ndark;
    int Nstar;
    double ri;
    double ro;
    double Mtot;
    double Mgas;
    double Mdark;
    double Mstar;
    } PA;

void usage(void);

int main(int argc, char **argv) {

    int i, j;
    int Nbin;
    int positionprecision;
    int gridtype;
    double r, dr, rmin, rmax, vol;
    double dx, dy, dz;
    double rcentre[3] = {0,0,0};
    TIPSY_HEADER th;
    GAS_PARTICLE gp;
    DARK_PARTICLE dp;
    STAR_PARTICLE sp;
    GAS_PARTICLE_DPP gpdpp;
    DARK_PARTICLE_DPP dpdpp;
    STAR_PARTICLE_DPP spdpp;
    PA *pa;
    XDR xdrs;

    positionprecision = 0; /* spp 0, dpp 1 */
    gridtype = 1; /* lin 0, log 1 */
    rmin = 0.1;
    rmax = 1;
    Nbin = 1;
    dr = 0;
    /*
    ** Read in arguments
    */
    i = 1;
    while (i < argc) {
	if (strcmp(argv[i],"-spp") == 0) {
            positionprecision = 0;
            i++;
            }
        else if (strcmp(argv[i],"-dpp") == 0) {
            positionprecision = 1;
            i++;
            }
        else if (strcmp(argv[i],"-lin") == 0) {
	    gridtype = 0;
            i++;
            }
        else if (strcmp(argv[i],"-log") == 0) {
	    gridtype = 1;
            i++;
            }
	else if (strcmp(argv[i],"-rmin") == 0) {
	    i++;
	    if (i >= argc) {
		usage();
		}
	    rmin= atof(argv[i]);
	    i++;
	    }
	else if (strcmp(argv[i],"-rmax") == 0) {
	    i++;
	    if (i >= argc) {
		usage();
		}
	    rmax = atof(argv[i]);
	    i++;
	    }
	else if (strcmp(argv[i],"-Nbin") == 0) {
	    i++;
	    if (i >= argc) {
		usage();
		}
	    Nbin = atof(argv[i]);
	    i++;
	    }
        else if (strcmp(argv[i],"-rxcen") == 0) {
            i++;
            if (i >= argc) {
                usage();
                }
            rcentre[0] = atof(argv[i]);
            i++;
            }
        else if (strcmp(argv[i],"-rycen") == 0) {
            i++;
            if (i >= argc) {
                usage();
                }
            rcentre[1] = atof(argv[i]);
            i++;
            }
        else if (strcmp(argv[i],"-rzcen") == 0) {
            i++;
            if (i >= argc) {
                usage();
                }
            rcentre[2] = atof(argv[i]);
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
    ** Initialise array
    */
    if (gridtype == 0) {
	dr = (rmax-rmin)/Nbin;
	}
    else if (gridtype == 1) {
	dr = (log(rmax)-log(rmin))/Nbin;
	}
    pa = malloc(Nbin*sizeof(PA));
    for (j = 0; j < Nbin; j++) {
	if (gridtype == 0) {
	    pa[j].ri = rmin + j*dr;
	    pa[j].ro = rmin + (j+1)*dr;
	    }
	else if (gridtype == 1) {
	    pa[j].ri = exp(log(rmin) + j*dr);
	    pa[j].ro = exp(log(rmin) + (j+1)*dr);
	    }
	pa[j].Ntot = 0;
	pa[j].Ngas = 0;
	pa[j].Ndark = 0;
	pa[j].Nstar = 0;
	pa[j].Mtot = 0;
	pa[j].Mgas = 0;
	pa[j].Mdark = 0;
	pa[j].Mstar = 0;
	}
    /*
    ** Read in particles and calculate profile
    */
    xdrstdio_create(&xdrs,stdin,XDR_DECODE);
    read_tipsy_standard_header(&xdrs,&th);
    if (positionprecision == 0) {
	for (i = 0; i < th.ngas; i++) {
	    read_tipsy_standard_gas(&xdrs,&gp);
	    }
	for (i = 0; i < th.ndark; i++) {
	    read_tipsy_standard_dark(&xdrs,&dp);
	    dx = dp.pos[0]-rcentre[0];
	    dy = dp.pos[1]-rcentre[1];
	    dz = dp.pos[2]-rcentre[2];
	    r = sqrt(dx*dy + dy*dy + dz*dz);
	    for (j = 0; j < Nbin; j++) {
		if ((pa[j].ri <= r) && (pa[j].ro > r)) {
		    pa[j].Ntot++;
		    pa[j].Ndark++;
		    pa[j].Mtot += dp.mass;
		    pa[j].Mdark += dp.mass;
		    }
		}
	    }
	for (i = 0; i < th.nstar; i++) {
	    read_tipsy_standard_star(&xdrs,&sp);
	    }
	}
    else if (positionprecision == 0) {
	for (i = 0; i < th.ngas; i++) {
	    read_tipsy_standard_gas_dpp(&xdrs,&gpdpp);
	    }
	for (i = 0; i < th.ndark; i++) {
	    read_tipsy_standard_dark_dpp(&xdrs,&dpdpp);
	    dx = dpdpp.pos[0]-rcentre[0];
	    dy = dpdpp.pos[1]-rcentre[1];
	    dz = dpdpp.pos[2]-rcentre[2];
	    r = sqrt(dx*dy + dy*dy + dz*dz);
	    for (j = 0; j < Nbin; j++) {
		if ((pa[j].ri <= r) && (pa[j].ro > r)) {
		    pa[j].Ntot++;
		    pa[j].Ndark++;
		    pa[j].Mtot += dp.mass;
		    pa[j].Mdark += dp.mass;
		    }
		}
	    }
	for (i = 0; i < th.nstar; i++) {
	    read_tipsy_standard_star_dpp(&xdrs,&spdpp);
	    }
	}
    /*
    ** Write output
    */
    for (j = 0; j < Nbin; j++) {
	if (gridtype == 0) {
	    fprintf(stdout,"%.6e %.6e %.6e ",pa[j].ri,(pa[j].ro+pa[j].ri)/2.0,pa[j].ro);
	    }
	else if (gridtype == 1) {
	    fprintf(stdout,"%.6e %.6e %.6e ",pa[j].ri,exp(log(pa[j].ri)+dr/2.0),pa[j].ro);
	    }
	fprintf(stdout,"%.6e %.6e %.6e %.6e ",pa[j].Mtot,pa[j].Mgas,pa[j].Mdark,pa[j].Mstar);
	vol = 4*M_PI*(pa[j].ro*pa[j].ro*pa[j].ro - pa[j].ri*pa[j].ri*pa[j].ri)/3.0;
	fprintf(stdout,"%.6e %.6e %.6e %.6e ",pa[j].Mtot/vol,pa[j].Mgas/vol,pa[j].Mdark/vol,pa[j].Mstar/vol);
	fprintf(stdout,"%d %d %d %d\n",pa[j].Ntot,pa[j].Ngas,pa[j].Ndark,pa[j].Nstar);
	}
    exit(0);
    }

void usage(void) {

    fprintf(stderr,"\n");
    fprintf(stderr,"Program calculates the profile of the input file around r_cen.\n");
    fprintf(stderr,"\n");
    fprintf(stderr,"You can specify the following arguments:\n");
    fprintf(stderr,"\n");
    fprintf(stderr,"-spp           : set this flag if input and output files have single precision positions (default)\n");
    fprintf(stderr,"-dpp           : set this flag if input and output files have double precision positions\n");
    fprintf(stderr,"-lin           : set this flag for linear grid\n");
    fprintf(stderr,"-log           : set this flag for logarithmic grid (default)\n");
    fprintf(stderr,"-rmin <value>  : minimum grid radius [LU]\n");
    fprintf(stderr,"-rmax <value>  : maximum grid radius [LU]\n");
    fprintf(stderr,"-Nbin <value>  : number of bins\n");
    fprintf(stderr,"-rxcen <value> : x-coordinate of centre [LU] (default: 0 LU)\n");
    fprintf(stderr,"-rycen <value> : y-coordinate of centre [LU] (default: 0 LU)\n");
    fprintf(stderr,"-rzcen <value> : z-coordinate of centre [LU] (default: 0 LU)\n");
    fprintf(stderr,"< <name>       : name of input file in tipsy standard binary format\n");
    fprintf(stderr,"> <name>       : name of output file\n");
    fprintf(stderr,"\n");
    exit(1);
    }
