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
    double veltot[3];
    double vel2tot[3];
    double velgas[3];
    double vel2gas[3];
    double veldark[3];
    double vel2dark[3];
    double velstar[3];
    double vel2star[3];
    } PA;

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

void usage(void);

int main(int argc, char **argv) {

    int i, j, k;
    int Nbin;
    int positionprecision;
    int gridtype;
    double r, dr, rmin, rmax, vol;
    double Menctot, Mencgas, Mencdark, Mencstar;
    int Nenctot, Nencgas, Nencdark, Nencstar;
    double pos[3], vel[3];
    double velproj[3], velmean, vel2mean;
    double erad[3], ephi[3], etheta[3];
    double rcentre[3] = {0,0,0}, vcentre[3] = {0,0,0};
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
    rmin = 0;
    rmax = 0;
    Nbin = 0;
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
        else if (strcmp(argv[i],"-vxcen") == 0) {
            i++;
            if (i >= argc) {
                usage();
                }
            vcentre[0] = atof(argv[i]);
            i++;
            }
        else if (strcmp(argv[i],"-vycen") == 0) {
            i++;
            if (i >= argc) {
                usage();
                }
            vcentre[1] = atof(argv[i]);
            i++;
            }
        else if (strcmp(argv[i],"-vzcen") == 0) {
            i++;
            if (i >= argc) {
                usage();
                }
            vcentre[2] = atof(argv[i]);
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
    ** Check some things
    */
    assert(Nbin > 0);
    assert(rmax > 0);
    assert(rmin >= 0);
    /*
    ** Initialise array
    */
    if (gridtype == 0) {
	dr = (rmax-rmin)/Nbin;
	}
    else if (gridtype == 1) {
	dr = (log(rmax)-log(rmin))/Nbin;
	}
    pa = malloc((Nbin+1)*sizeof(PA));
    for (j = 0; j < (Nbin+1); j++) {
	if (gridtype == 0) {
	    if (j == 0) {
		pa[j].ri = 0;
		pa[j].ro = rmin;
		}
	    else {
		pa[j].ri = rmin + (j-1)*dr;
		pa[j].ro = rmin + j*dr;
		}
	    }
	else if (gridtype == 1) {
	    if (j == 0) {
		pa[j].ri = 0;
		pa[j].ro = rmin;
		}
	    else {
		pa[j].ri = exp(log(rmin) + (j-1)*dr);
		pa[j].ro = exp(log(rmin) + j*dr);
		}
	    }
	pa[j].Ntot = 0;
	pa[j].Ngas = 0;
	pa[j].Ndark = 0;
	pa[j].Nstar = 0;
	pa[j].Mtot = 0;
	pa[j].Mgas = 0;
	pa[j].Mdark = 0;
	pa[j].Mstar = 0;
	for (i = 0; i < 3; i++) {
	    pa[j].veltot[i] = 0;
	    pa[j].velgas[i] = 0;
	    pa[j].veldark[i] = 0;
	    pa[j].velstar[i] = 0;
	    pa[j].vel2tot[i] = 0;
	    pa[j].vel2gas[i] = 0;
	    pa[j].vel2dark[i] = 0;
	    pa[j].vel2star[i] = 0;
	    }
	}
    /*
    ** Read in particles and calculate profile
    */
    xdrstdio_create(&xdrs,stdin,XDR_DECODE);
    read_tipsy_standard_header(&xdrs,&th);
    if (positionprecision == 0) {
	for (i = 0; i < th.ngas; i++) {
	    read_tipsy_standard_gas(&xdrs,&gp);
	    for (j = 0; j < 3; j++) {
		pos[j] = gp.pos[j]-rcentre[j];
		vel[j] = gp.vel[j]-vcentre[j];
		}
	    r = sqrt(pos[0]*pos[0]+pos[1]*pos[1]+pos[2]*pos[2]);
	    for (j = 0; j < (Nbin+1); j++) {
		if ((pa[j].ri <= r) && (pa[j].ro > r)) {
		    calculate_unit_vectors(pos,erad,ephi,etheta);
		    velproj[0] = vel[0]*erad[0]+vel[1]*erad[1]+vel[2]*erad[2];
		    velproj[1] = vel[0]*ephi[0]+vel[1]*ephi[1]+vel[2]*ephi[2];
		    velproj[2] = vel[0]*etheta[0]+vel[1]*etheta[1]+vel[2]*etheta[2];
		    pa[j].Ntot++;
		    pa[j].Ngas++;
		    pa[j].Mtot += gp.mass;
		    pa[j].Mgas += gp.mass;
		    for (k = 0; k < 3; k++) {
			pa[j].veltot[k] += velproj[k];
			pa[j].velgas[k] += velproj[k];
			pa[j].vel2tot[k] += velproj[k]*velproj[k];
			pa[j].vel2gas[k] += velproj[k]*velproj[k];
			}
		    }
		}
	    }
	for (i = 0; i < th.ndark; i++) {
	    read_tipsy_standard_dark(&xdrs,&dp);
	    for (j = 0; j < 3; j++) {
		pos[j] = dp.pos[j]-rcentre[j];
		vel[j] = dp.vel[j]-vcentre[j];
		}
	    r = sqrt(pos[0]*pos[0]+pos[1]*pos[1]+pos[2]*pos[2]);
	    for (j = 0; j < (Nbin+1); j++) {
		if ((pa[j].ri <= r) && (pa[j].ro > r)) {
		    calculate_unit_vectors(pos,erad,ephi,etheta);
		    velproj[0] = vel[0]*erad[0]+vel[1]*erad[1]+vel[2]*erad[2];
		    velproj[1] = vel[0]*ephi[0]+vel[1]*ephi[1]+vel[2]*ephi[2];
		    velproj[2] = vel[0]*etheta[0]+vel[1]*etheta[1]+vel[2]*etheta[2];
		    pa[j].Ntot++;
		    pa[j].Ndark++;
		    pa[j].Mtot += dp.mass;
		    pa[j].Mdark += dp.mass;
		    for (k = 0; k < 3; k++) {
			pa[j].veltot[k] += velproj[k];
			pa[j].veldark[k] += velproj[k];
			pa[j].vel2tot[k] += velproj[k]*velproj[k];
			pa[j].vel2dark[k] += velproj[k]*velproj[k];
			}
		    }
		}
	    }
	for (i = 0; i < th.nstar; i++) {
	    read_tipsy_standard_star(&xdrs,&sp);
	    read_tipsy_standard_dark(&xdrs,&dp);
	    for (j = 0; j < 3; j++) {
		pos[j] = sp.pos[j]-rcentre[j];
		vel[j] = sp.vel[j]-vcentre[j];
		}
	    r = sqrt(pos[0]*pos[0]+pos[1]*pos[1]+pos[2]*pos[2]);
	    for (j = 0; j < (Nbin+1); j++) {
		if ((pa[j].ri <= r) && (pa[j].ro > r)) {
		    calculate_unit_vectors(pos,erad,ephi,etheta);
		    velproj[0] = vel[0]*erad[0]+vel[1]*erad[1]+vel[2]*erad[2];
		    velproj[1] = vel[0]*ephi[0]+vel[1]*ephi[1]+vel[2]*ephi[2];
		    velproj[2] = vel[0]*etheta[0]+vel[1]*etheta[1]+vel[2]*etheta[2];
		    pa[j].Ntot++;
		    pa[j].Nstar++;
		    pa[j].Mtot += sp.mass;
		    pa[j].Mstar += sp.mass;
		    for (k = 0; k < 3; k++) {
			pa[j].veltot[k] += velproj[k];
			pa[j].velstar[k] += velproj[k];
			pa[j].vel2tot[k] += velproj[k]*velproj[k];
			pa[j].vel2star[k] += velproj[k]*velproj[k];
			}
		    }
		}
	    }
	}
    else if (positionprecision == 1) {
	for (i = 0; i < th.ngas; i++) {
	    read_tipsy_standard_gas_dpp(&xdrs,&gpdpp);
	    for (j = 0; j < 3; j++) {
		pos[j] = gpdpp.pos[j]-rcentre[j];
		vel[j] = gpdpp.vel[j]-vcentre[j];
		}
	    r = sqrt(pos[0]*pos[0]+pos[1]*pos[1]+pos[2]*pos[2]);
	    for (j = 0; j < (Nbin+1); j++) {
		if ((pa[j].ri <= r) && (pa[j].ro > r)) {
		    calculate_unit_vectors(pos,erad,ephi,etheta);
		    velproj[0] = vel[0]*erad[0]+vel[1]*erad[1]+vel[2]*erad[2];
		    velproj[1] = vel[0]*ephi[0]+vel[1]*ephi[1]+vel[2]*ephi[2];
		    velproj[2] = vel[0]*etheta[0]+vel[1]*etheta[1]+vel[2]*etheta[2];
		    pa[j].Ntot++;
		    pa[j].Ngas++;
		    pa[j].Mtot += gpdpp.mass;
		    pa[j].Mgas += gpdpp.mass;
		    for (k = 0; k < 3; k++) {
			pa[j].veltot[k] += velproj[k];
			pa[j].velgas[k] += velproj[k];
			pa[j].vel2tot[k] += velproj[k]*velproj[k];
			pa[j].vel2gas[k] += velproj[k]*velproj[k];
			}
		    }
		}
	    }
	for (i = 0; i < th.ndark; i++) {
	    read_tipsy_standard_dark_dpp(&xdrs,&dpdpp);
	    for (j = 0; j < 3; j++) {
		pos[j] = dpdpp.pos[j]-rcentre[j];
		vel[j] = dpdpp.vel[j]-vcentre[j];
		}
	    r = sqrt(pos[0]*pos[0]+pos[1]*pos[1]+pos[2]*pos[2]);
	    for (j = 0; j < (Nbin+1); j++) {
		if ((pa[j].ri <= r) && (pa[j].ro > r)) {
		    calculate_unit_vectors(pos,erad,ephi,etheta);
		    velproj[0] = vel[0]*erad[0]+vel[1]*erad[1]+vel[2]*erad[2];
		    velproj[1] = vel[0]*ephi[0]+vel[1]*ephi[1]+vel[2]*ephi[2];
		    velproj[2] = vel[0]*etheta[0]+vel[1]*etheta[1]+vel[2]*etheta[2];
		    pa[j].Ntot++;
		    pa[j].Ndark++;
		    pa[j].Mtot += dpdpp.mass;
		    pa[j].Mdark += dpdpp.mass;
		    for (k = 0; k < 3; k++) {
			pa[j].veltot[k] += velproj[k];
			pa[j].veldark[k] += velproj[k];
			pa[j].vel2tot[k] += velproj[k]*velproj[k];
			pa[j].vel2dark[k] += velproj[k]*velproj[k];
			}
		    }
		}
	    }
	for (i = 0; i < th.nstar; i++) {
	    read_tipsy_standard_star_dpp(&xdrs,&spdpp);
	    for (j = 0; j < 3; j++) {
		pos[j] = spdpp.pos[j]-rcentre[j];
		vel[j] = spdpp.vel[j]-vcentre[j];
		}
	    r = sqrt(pos[0]*pos[0]+pos[1]*pos[1]+pos[2]*pos[2]);
	    for (j = 0; j < (Nbin+1); j++) {
		if ((pa[j].ri <= r) && (pa[j].ro > r)) {
		    calculate_unit_vectors(pos,erad,ephi,etheta);
		    velproj[0] = vel[0]*erad[0]+vel[1]*erad[1]+vel[2]*erad[2];
		    velproj[1] = vel[0]*ephi[0]+vel[1]*ephi[1]+vel[2]*ephi[2];
		    velproj[2] = vel[0]*etheta[0]+vel[1]*etheta[1]+vel[2]*etheta[2];
		    pa[j].Ntot++;
		    pa[j].Nstar++;
		    pa[j].Mtot += spdpp.mass;
		    pa[j].Mstar += spdpp.mass;
		    for (k = 0; k < 3; k++) {
			pa[j].veltot[k] += velproj[k];
			pa[j].velstar[k] += velproj[k];
			pa[j].vel2tot[k] += velproj[k]*velproj[k];
			pa[j].vel2star[k] += velproj[k]*velproj[k];
			}
		    }
		}
	    }
	}
    /*
    ** Write output
    */
    for (j = 0; j < (Nbin+1); j++) {
	vol = 4*M_PI*(pa[j].ro*pa[j].ro*pa[j].ro - pa[j].ri*pa[j].ri*pa[j].ri)/3.0;
	if (j == 0) {
	    fprintf(stdout,"%.6e %.6e %.6e %.6e ",pa[j].ri,(pa[j].ro+pa[j].ri)/2.0,pa[j].ro,vol);
	    }
	else {
	    if (gridtype == 0) {
		fprintf(stdout,"%.6e %.6e %.6e %.6e ",pa[j].ri,(pa[j].ro+pa[j].ri)/2.0,pa[j].ro,vol);
		}
	    else if (gridtype == 1) {
		fprintf(stdout,"%.6e %.6e %.6e %.6e ",pa[j].ri,exp(log(pa[j].ri)+dr/2.0),pa[j].ro,vol);
		}
	    }
	fprintf(stdout,"%.6e %.6e %.6e %.6e ",pa[j].Mtot,pa[j].Mgas,pa[j].Mdark,pa[j].Mstar);
	Menctot = 0;
	Mencgas = 0;
	Mencdark = 0;
	Mencstar = 0;
	for (i = 0; i <= j; i++) {
	    Menctot += pa[i].Mtot;
	    Mencgas += pa[i].Mgas;
	    Mencdark += pa[i].Mdark;
	    Mencstar += pa[i].Mstar;
	    }
	fprintf(stdout,"%.6e %.6e %.6e %.6e ",Menctot,Mencgas,Mencdark,Mencstar);
	fprintf(stdout,"%.6e %.6e %.6e %.6e ",pa[j].Mtot/vol,pa[j].Mgas/vol,pa[j].Mdark/vol,pa[j].Mstar/vol);
	fprintf(stdout,"%d %d %d %d ",pa[j].Ntot,pa[j].Ngas,pa[j].Ndark,pa[j].Nstar);
	Nenctot = 0;
	Nencgas = 0;
	Nencdark = 0;
	Nencstar = 0;
	for (i = 0; i <= j; i++) {
	    Nenctot += pa[i].Ntot;
	    Nencgas += pa[i].Ngas;
	    Nencdark += pa[i].Ndark;
	    Nencstar += pa[i].Nstar;
	    }
	fprintf(stdout,"%d %d %d %d ",Nenctot,Nencgas,Nencdark,Nencstar);
	fprintf(stdout,"%.6e %.6e %.6e %.6e ",pa[j].Ntot/vol,pa[j].Ngas/vol,pa[j].Ndark/vol,pa[j].Nstar/vol);
	for (i = 0; i < 3; i++) {
	    velmean = pa[j].veltot[i]/pa[j].Ntot;
	    vel2mean = pa[j].vel2tot[i]/pa[j].Ntot;
	    fprintf(stdout,"%.6e %.6e ",velmean,sqrt(vel2mean-velmean*velmean));
	    }
	for (i = 0; i < 3; i++) {
	    velmean = pa[j].velgas[i]/pa[j].Ngas;
	    vel2mean = pa[j].vel2gas[i]/pa[j].Ngas;
	    fprintf(stdout,"%.6e %.6e ",velmean,sqrt(vel2mean-velmean*velmean));
	    }
	for (i = 0; i < 3; i++) {
	    velmean = pa[j].veldark[i]/pa[j].Ndark;
	    vel2mean = pa[j].vel2dark[i]/pa[j].Ndark;
	    fprintf(stdout,"%.6e %.6e ",velmean,sqrt(vel2mean-velmean*velmean));
	    }
	for (i = 0; i < 3; i++) {
	    velmean = pa[j].velstar[i]/pa[j].Nstar;
	    vel2mean = pa[j].vel2star[i]/pa[j].Nstar;
	    fprintf(stdout,"%.6e %.6e ",velmean,sqrt(vel2mean-velmean*velmean));
	    }
	fprintf(stdout,"\n");
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
    fprintf(stderr,"-Nbin <value>  : number of bins between rmin and rmax\n");
    fprintf(stderr,"-rxcen <value> : x-coordinate of centre [LU] (default: 0 LU)\n");
    fprintf(stderr,"-rycen <value> : y-coordinate of centre [LU] (default: 0 LU)\n");
    fprintf(stderr,"-rzcen <value> : z-coordinate of centre [LU] (default: 0 LU)\n");
    fprintf(stderr,"< <name>       : name of input file in tipsy standard binary format\n");
    fprintf(stderr,"> <name>       : name of output file\n");
    fprintf(stderr,"\n");
    exit(1);
    }
