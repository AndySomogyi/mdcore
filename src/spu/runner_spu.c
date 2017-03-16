/*******************************************************************************
 * This file is part of mdcore.
 * Coypright (c) 2010 Pedro Gonnet (gonnet@maths.ox.ac.uk)
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 * 
 ******************************************************************************/

/* make sure we're in cell mode (so the headers won't get messed-up). */
#ifndef CELL
    #define CELL
#endif

/* include some header files */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <spu_mfcio.h>
#include <spu_intrinsics.h>
#include <libmisc.h>
#include <simdmath.h>
#include "load_vec_float4.h"
#include "dot_product3.h"

/* some definitions */
#define bitesize                            3
#define maxparts                            256
#define maxqstack                           128
#define vec_splat(_a, _b)	                spu_splats(spu_extract(_a, _b))

/* potential flags */
#define potential_flag_none                 0
#define potential_flag_LJ126                1
#define potential_flag_Ewald                2
#define potential_flag_Coulomb              4

/* potential degree (this has to match potential.h!) */
#define potential_degree                    5
#define potential_chunk                     (potential_degree+3)

/* particle flags */
#define part_flag_none                      0
#define part_flag_frozen                    1
#define part_flag_ghost                     2

/* declare the LOCAL types we will use */
struct potential {
    int n;
    unsigned int flags;
    float alpha[3];
    float data[];
    };
    
/* This has to match part.h! note the use of vector types instead
   of arrays of FPTYPE... */
struct part {
    vector float x;
    vector float v;
    vector float f;
    float q;
    int id, vid;
    short int type;
    unsigned short int flags;
    };
    
struct data {
    unsigned long long cell_addr;
    int nr_cells;
    float cutoff;
    float h[3];
    int max_type;
    struct potential *p[];
    };
    
struct cell {
    int ni;
    unsigned long long ai;
    };
    
struct data *data = NULL;
unsigned int rcount = 0;


/*////////////////////////////////////////////////////////////////////////////// */
/* void potential_eval_vec */
//
/* evaluates the given potential at the given point. */
/*////////////////////////////////////////////////////////////////////////////// */

inline void potential_eval_vec ( struct potential *p[4] , vector float r2 , vector float *e , vector float *f ) {

    int i, k;
    vector unsigned int ind;
    vector float x, ee, eff, r, alpha0, alpha1, alpha2, c, mi, hi;
    float *data[4];
    
    /* get the sqrt of r2 */
    r = sqrtf4( r2 );
    
    /* load the alphas */
    alpha0 = _load_vec_float4( p[0]->alpha[0] , p[1]->alpha[0] , p[2]->alpha[0] , p[3]->alpha[0] );
    alpha1 = _load_vec_float4( p[0]->alpha[1] , p[1]->alpha[1] , p[2]->alpha[1] , p[3]->alpha[1] );
    alpha2 = _load_vec_float4( p[0]->alpha[2] , p[1]->alpha[2] , p[2]->alpha[2] , p[3]->alpha[2] );
    
    /* compute the index (spu_convtu returns zero for negative values). */
    ind = spu_convtu( spu_madd( r , spu_madd( r , alpha2 , alpha1 ) , alpha0 ) , 0 );
    
    /* get a pointer to the data for this interval */
    for ( k = 0 ; k < 4 ; k++ )
        data[k] = &( p[k]->data[ potential_chunk * spu_extract( ind , k ) ] );
    
    /* get mi and hi */
    mi = _load_vec_float4( data[0][0] , data[1][0] , data[2][0] , data[3][0] );
    hi = _load_vec_float4( data[0][1] , data[1][1] , data[2][1] , data[3][1] );
        
    /* adjust x to the interval */
    x = spu_mul( spu_sub( r , mi ) , hi );
    
    /* compute the potential and its derivative */
    eff = _load_vec_float4( data[0][2] , data[1][2] , data[2][2] , data[3][2] );
    ee = spu_mul( eff , x );
    c = _load_vec_float4( data[0][3] , data[1][3] , data[2][3] , data[3][3] );
    ee = spu_add( ee , c );
    for ( i = 4 ; i < potential_degree+3 ; i++ ) {
        eff = spu_madd( eff , x , ee );
        c = _load_vec_float4( data[0][i] , data[1][i] , data[2][i] , data[3][i] );
        ee = spu_madd( ee , x , c );
        }

    /* store the result */
    *e = ee; *f = eff * ( hi / r );    
    
    }


/*////////////////////////////////////////////////////////////////////////////// */
/* void potential_eval_vec2 */
//
/* evaluates the given potentials at the given point. */
/*////////////////////////////////////////////////////////////////////////////// */

inline void potential_eval_vec2 ( struct potential *p[8] , vector float *r2 , vector float *e , vector float *f ) {

    int i, k;
    vector unsigned int ind_1, ind_2;
    vector float x_1, ee_1, eff_1, r_1, alpha0_1, alpha1_1, alpha2_1, c_1, mi_1, hi_1;
    vector float x_2, ee_2, eff_2, r_2, alpha0_2, alpha1_2, alpha2_2, c_2, mi_2, hi_2;
    float *data[8];
    
    /* get the sqrt of r2 */
    r_1 = sqrtf4( r2[0] );
    r_2 = sqrtf4( r2[1] );
    
    /* load the alphas */
    alpha0_1 = _load_vec_float4( p[0]->alpha[0] , p[1]->alpha[0] , p[2]->alpha[0] , p[3]->alpha[0] );
    alpha1_1 = _load_vec_float4( p[0]->alpha[1] , p[1]->alpha[1] , p[2]->alpha[1] , p[3]->alpha[1] );
    alpha2_1 = _load_vec_float4( p[0]->alpha[2] , p[1]->alpha[2] , p[2]->alpha[2] , p[3]->alpha[2] );
    alpha0_2 = _load_vec_float4( p[4]->alpha[0] , p[5]->alpha[0] , p[6]->alpha[0] , p[7]->alpha[0] );
    alpha1_2 = _load_vec_float4( p[4]->alpha[1] , p[5]->alpha[1] , p[6]->alpha[1] , p[7]->alpha[1] );
    alpha2_2 = _load_vec_float4( p[4]->alpha[2] , p[5]->alpha[2] , p[6]->alpha[2] , p[7]->alpha[2] );
    
    /* compute the index */
    ind_1 = spu_convtu( spu_madd( r_1 , spu_madd( r_1 , alpha2_1 , alpha1_1 ) , alpha0_1 ) , 0 );
    ind_2 = spu_convtu( spu_madd( r_2 , spu_madd( r_2 , alpha2_2 , alpha1_2 ) , alpha0_2 ) , 0 );
    
    /* get a pointer to the data for this interval */
    for ( k = 0 ; k < 4 ; k++ )
        data[k] = &( p[k]->data[ potential_chunk * spu_extract( ind_1 , k ) ] );
    for ( k = 0 ; k < 4 ; k++ )
        data[k+4] = &( p[k+4]->data[ potential_chunk * spu_extract( ind_2 , k ) ] );
    
    /* get mi and hi */
    mi_1 = _load_vec_float4( data[0][0] , data[1][0] , data[2][0] , data[3][0] );
    hi_1 = _load_vec_float4( data[0][1] , data[1][1] , data[2][1] , data[3][1] );
    mi_2 = _load_vec_float4( data[4][0] , data[5][0] , data[6][0] , data[7][0] );
    hi_2 = _load_vec_float4( data[4][1] , data[5][1] , data[6][1] , data[7][1] );
        
    /* adjust x to the interval */
    x_1 = spu_mul( spu_sub( r_1 , mi_1 ) , hi_1 );
    x_2 = spu_mul( spu_sub( r_2 , mi_2 ) , hi_2 );
    
    /* compute the potential and its derivative */
    eff_1 = _load_vec_float4( data[0][2] , data[1][2] , data[2][2] , data[3][2] );
    eff_2 = _load_vec_float4( data[4][2] , data[5][2] , data[6][2] , data[7][2] );
    ee_1 = spu_mul( eff_1 , x_1 );
    ee_2 = spu_mul( eff_2 , x_2 );
    c_1 = _load_vec_float4( data[0][3] , data[1][3] , data[2][3] , data[3][3] );
    c_2 = _load_vec_float4( data[4][3] , data[5][3] , data[6][3] , data[7][3] );
    ee_1 = spu_add( ee_1 , c_1 );
    ee_2 = spu_add( ee_2 , c_2 );
    for ( i = 4 ; i < potential_degree+3 ; i++ ) {
        eff_1 = spu_madd( eff_1 , x_1 , ee_1 );
        eff_2 = spu_madd( eff_2 , x_2 , ee_2 );
        c_1 = _load_vec_float4( data[0][i] , data[1][i] , data[2][i] , data[3][i] );
        c_2 = _load_vec_float4( data[4][i] , data[5][i] , data[6][i] , data[7][i] );
        ee_1 = spu_madd( ee_1 , x_1 , c_1 );
        ee_2 = spu_madd( ee_2 , x_2 , c_2 );
        }

    /* store the result */
    e[0] = ee_1; f[0] = eff_1 * ( hi_1 / r_1 );    
    e[1] = ee_2; f[1] = eff_2 * ( hi_2 / r_2 );    
    
    }


/*////////////////////////////////////////////////////////////////////////////// */
/* void potential_eval_expl */
//
/* evaluates the given potential at the given point. */
/*////////////////////////////////////////////////////////////////////////////// */

#ifdef USE_DOUBLES
inline void potential_eval_expl ( struct potential *p , double r2 , double *e , double *f ) {

    const double isqrtpi = 0.56418958354775628695;
    const double kappa = 3.0;
    double r = sqrt(r2), ir = 1.0 / r, ir2 = ir * ir, ir4, ir6, ir12, t1, t2;

    /* init f and e */
    *e = 0.0; *f = 0.0;

    /* do we have a Lennard-Jones interaction? */
    if ( p->flags & potential_flag_LJ126 ) {
    
        /* init some variables */
        ir4 = ir2 * ir2; ir6 = ir4 * ir2; ir12 = ir6 * ir6;
        
        /* compute the energy and the force */
        *e = ( p->alpha[0] * ir12 - p->alpha[1] * ir6 );
        *f = 6.0 * ir * ( -2.0 * p->alpha[0] * ir12 + p->alpha[1] * ir6 );
    
        }
        
    /* do we have an Ewald short-range part? */
    if ( p->flags & potential_flag_Ewald ) {
    
        /* get some values we will re-use */
        t2 = r * kappa;
        t1 = erfc( t2 );
    
        /* compute the energy and the force */
        *e += p->alpha[2] * t1 * ir;
        *f += p->alpha[2] * ( -2.0 * isqrtpi * exp( -t2 * t2 ) * kappa * ir - t1 * ir2 );
    
        }
    
    /* do we have a Coulomb interaction? */
    if ( p->flags & potential_flag_Coulomb ) {
    
        /* get some values we will re-use */
        t2 = r * kappa;
        t1 = erfc( t2 );
    
        /* compute the energy and the force */
        *e += p->alpha[2] * ir;
        *f += -p->alpha[2] * ir2;
    
        }
    
    }
#else
inline void potential_eval_expl ( struct potential *p , float r2 , float *e , float *f ) {

    const float isqrtpi = 0.56418958354775628695f;
    const float kappa = 3.0f;
    float r = sqrt(r2), ir = 1.0f / r, ir2 = ir * ir, ir4, ir6, ir12, t1, t2;

    /* init f and e */
    *e = 0.0f; *f = 0.0f;

    /* do we have a Lennard-Jones interaction? */
    if ( p->flags & potential_flag_LJ126 ) {
    
        /* init some variables */
        ir4 = ir2 * ir2; ir6 = ir4 * ir2; ir12 = ir6 * ir6;
        
        /* compute the energy and the force */
        *e = ( p->alpha[0] * ir12 - p->alpha[1] * ir6 );
        *f = 6.0f * ir * ( -2.0f * p->alpha[0] * ir12 + p->alpha[1] * ir6 );
    
        }
        
    /* do we have an Ewald short-range part? */
    if ( p->flags & potential_flag_Ewald ) {
    
        /* get some values we will re-use */
        t2 = r * kappa;
        t1 = erfc( t2 );
    
        /* compute the energy and the force */
        *e += p->alpha[2] * t1 * ir;
        *f += p->alpha[2] * ( -2.0f * isqrtpi * exp( -t2 * t2 ) * kappa * ir - t1 * ir2 );
    
        }
    
    /* do we have a Coulomb interaction? */
    if ( p->flags & potential_flag_Coulomb ) {
    
        /* get some values we will re-use */
        t2 = r * kappa;
        t1 = erfc( t2 );
    
        /* compute the energy and the force */
        *e += p->alpha[2] * ir;
        *f += -p->alpha[2] * ir2;
    
        }
    
    }
#endif


/*////////////////////////////////////////////////////////////////////////////// */
/* void potential_eval */
//
/* evaluates the given potential at the given point. */
/*////////////////////////////////////////////////////////////////////////////// */

#ifdef USE_DOUBLES
inline void potential_eval ( struct potential *p , double r2 , double *e , double *f ) {

    int ind, k;
    double x, ee, eff, r = sqrt(r2);
    float *data;
    
    /* is r in the house? */
    /* if ( r < p->a || r > p->b ) */
    /*     printf("potential_eval: requested potential at r=%e, not in [%e,%e].\n",r,p->a,p->b); */
    
    /* compute the index */
    ind = fmax( 0.0 , p->alpha[0] + r * (p->alpha[1] + r * p->alpha[2]) );
    
    /* get a pointer to the data for this interval */
    data = &( p->data[ potential_chunk * ind ] );
    
    /* adjust x to the interval */
    x = (r - data[0]) * data[1];
    
    /* compute the potential and its derivative */
    ee = data[2] * x + data[3];
    eff = data[2];
    for ( k = 4 ; k < potential_chunk ; k++ ) {
        eff = eff * x + ee;
        ee = ee * x + data[k];
        }

    /* store the result */
    *e = ee; *f = eff * data[1] / r;
    
    }
#else
inline void potential_eval ( struct potential *p , float r2 , float *e , float *f ) {

    int ind, k;
    float x, ee, eff, *data, r = sqrt(r2);
    
    /* compute the index */
    ind = fmaxf( 0.0f , p->alpha[0] + r * (p->alpha[1] + r * p->alpha[2]) );
    
    /* get a pointer to the data for this interval */
    data = &( p->data[ potential_chunk * ind ] );
    
    /* adjust x to the interval */
    x = (r - data[0]) * data[1];
    
    /* compute the potential and its derivative */
    ee = data[2] * x + data[3];
    eff = data[2];
    for ( k = 4 ; k < potential_chunk ; k++ ) {
        eff = eff * x + ee;
        ee = ee * x + data[k];
        }

    /* store the result */
    *e = ee; *f = eff * data[1] / r;
    
    }
#endif


#ifdef VECTORIZE
/*////////////////////////////////////////////////////////////////////////////// */
/* sortedpair */
//
/* compute the pairwise interactions for the given pair using the sorted */
/* interactions algorithm */
/*////////////////////////////////////////////////////////////////////////////// */

void sortedpair ( int ni , int nj , struct part *pi , struct part *pj , float *shift ) {

    struct part *part_i, *part_j;
    struct potential *pot, *potv[8];
    float cutoff, cutoff2;
    float d[2*maxparts], temp, pivot;
    int ind[2*maxparts], left[2*maxparts], count = 0, lcount = 0, pcount = 0;
    int i, j, k, imax, qpos, lo, hi, pjoff, emt;
    struct {
        short int lo, hi;
        } qstack[maxqstack];
    float r2;
    #ifdef VEC2
        vector float *effi[8], *effj[8], dx, dxv[8], fv[2], ev[2], r2v[2], nshiftv, shiftv, tempv, pjx;
    #else
        vector float *effi[4], *effj[4], dx, dxv[4], fv, ev, r2v, nshiftv, shiftv, tempv, pjx;
    #endif
    
    /* get the space and cutoff */
    cutoff = data->cutoff;
    cutoff2 = cutoff * cutoff;
    emt = data->max_type;
        
    /* init r2v and make the compiler happy */
    #ifdef VEC2
        r2v[0] = spu_splats( 0.0f );
        r2v[1] = spu_splats( 0.0f );
    #else
        r2v = spu_splats( 0.0f );
    #endif
    
    /* extract the shift vector */
    shiftv = _load_vec_float4( shift[0] , shift[1] , shift[2] , 0.0f );
        
    /* start by filling the particle ids of both cells into ind and d */
    tempv = spu_mul( shiftv , shiftv );
    nshiftv = spu_mul( shiftv , rsqrtf4( spu_splats( spu_extract( tempv , 0 ) + spu_extract( tempv , 1 ) + spu_extract( tempv , 2 ) ) ) );
    for ( i = 0 ; i < ni ; i++ ) {
        part_i = &( pi[i] );
        ind[count] = -i - 1;
        tempv = spu_mul( nshiftv , part_i->x );
        d[count] = spu_extract( tempv , 0 ) + spu_extract( tempv , 1 ) + spu_extract( tempv , 2 );
        /* d[count] = _dot_product3( nshiftv , part_i->x ); */
        count += 1;
        }
    for ( i = 0 ; i < nj ; i++ ) {
        part_i = &( pj[i] );
        ind[count] = i;
        tempv = spu_mul( nshiftv , spu_add( part_i->x , shiftv ) );
        d[count] = spu_extract( tempv , 0 ) + spu_extract( tempv , 1 ) + spu_extract( tempv , 2 ) - cutoff;
        /* d[count] = _dot_product3( nshiftv , spu_add( part_i->x , shiftv ) ) - cutoff; */
        count += 1;
        }
        
    /* sort with quicksort */
    qstack[0].lo = 0; qstack[0].hi = count - 1; qpos = 0;
    while ( qpos >= 0 ) {
        lo = qstack[qpos].lo; hi = qstack[qpos].hi;
        qpos -= 1;
        if ( hi - lo < 10 ) {
            for ( i = lo ; i < hi ; i++ ) {
                imax = i;
                for ( j = i+1 ; j <= hi ; j++ )
                    if ( __builtin_expect( d[j] > d[imax] , 0 ) )
                        imax = j;
                if ( __builtin_expect( imax != i , 1 ) ) {
                    k = ind[imax]; ind[imax] = ind[i]; ind[i] = k;
                    temp = d[imax]; d[imax] = d[i]; d[i] = temp;
                    }
                }
            }
        else {
            pivot = d[ ( lo + hi ) / 2 ];
            i = lo; j = hi;
            while ( i <= j ) {
                while ( d[i] > pivot ) i++;
                while ( d[j] < pivot ) j--;
                if ( i <= j ) {
                    if ( __builtin_expect( i < j , 1 ) ) {
                        k = ind[i]; ind[i] = ind[j]; ind[j] = k;
                        temp = d[i]; d[i] = d[j]; d[j] = temp;
                        }
                    i += 1; j -= 1;
                    }
                }
            if ( __builtin_expect( lo < j , 1 ) ) {
                qpos += 1;
                qstack[qpos].lo = lo;
                qstack[qpos].hi = j;
                }
            if ( __builtin_expect( i < hi , 1 ) ) {
                qpos += 1;
                qstack[qpos].lo = i;
                qstack[qpos].hi = hi;
                }
            }
        }
        
    /* sort with selection sort */
    /* for ( i = 0 ; i < count-1 ; i++ ) { */
    /*     imax = i; */
    /*     for ( j = i+1 ; j < count ; j++ ) */
    /*         if ( __builtin_expect( d[j] > d[imax] , 0 ) ) */
    /*             imax = j; */
    /*     if ( __builtin_expect( imax != i , 1 ) ) { */
    /*         k = ind[imax]; ind[imax] = ind[i]; ind[i] = k; */
    /*         temp = d[imax]; d[imax] = d[i]; d[i] = temp; */
    /*         } */
    /*     } */
        
    /* loop over the sorted list of particles */
    for ( i = 0 ; i < count ; i++ ) {
    
        /* is this a particle from the left? */
        if ( ind[i] < 0 )
            left[lcount++] = -ind[i] - 1;
            
        /* it's from the right, interact with all left particles */
        else {
        
            /* get a handle on this particle */
            part_j = &( pj[ind[i]] );
            pjx = spu_add( part_j->x , shiftv );
            pjoff = emt * part_j->type;
        
            /* loop over the left particles */
            for ( j = lcount-1 ; j >= 0 ; j-- ) {
            
                /* get a handle on the second particle */
                part_i = &( pi[left[j]] );
            
                /* get the distance between both particles */
                dx = spu_sub( part_i->x , pjx );
                r2 = _dot_product3( dx , dx );
                
                /* is this within cutoff? */
                if ( r2 > cutoff2 )
                    continue;
                
                /* fetch the potential, if any */
                pot = data->p[ pjoff + part_i->type ];
                if ( __builtin_expect( pot == NULL , 0 ) )
                    continue;
                    
                /* add this interaction to the interaction queue */
                #ifdef VEC2
                    r2v[pcount/4][pcount%4] = r2;
                    dxv[pcount] = dx;
                    effi[pcount] = &( part_i->f );
                    effj[pcount] = &( part_j->f );
                    potv[pcount] = pot;
                    pcount += 1;
                    // rcount += 1;

                    /* do we have a full set to evaluate? */
                    if ( __builtin_expect( pcount == 8 , 0 ) ) {

                        /* evaluate the potentials */
                        potential_eval_vec2( potv , r2v , ev , fv );

                        /* for each entry, update the forces */
                        for ( k = 0 ; k < 8 ; k++ ) {
                            tempv = spu_mul( dxv[k] , vec_splat( fv[k/4] , k%4 ) );
                            *effi[k] = spu_sub( *effi[k] , tempv );
                            *effj[k] = spu_add( *effj[k] , tempv );
                            }

                        /* re-set the counter */
                        pcount = 0;

                        }
                #else
                    r2v[pcount] = r2;
                    dxv[pcount] = dx;
                    effi[pcount] = &( part_i->f );
                    effj[pcount] = &( part_j->f );
                    potv[pcount] = pot;
                    pcount += 1;
                    // rcount += 1;

                    /* do we have a full set to evaluate? */
                    if ( __builtin_expect( pcount == 4 , 0 ) ) {

                        /* evaluate the potentials */
                        potential_eval_vec( potv , r2v , &ev , &fv );

                        /* for each entry, update the forces */
                        for ( k = 0 ; k < 4 ; k++ ) {
                            tempv = spu_mul( dxv[k] , vec_splat( fv , k ) );
                            *effi[k] = spu_sub( *effi[k] , tempv );
                            *effj[k] = spu_add( *effj[k] , tempv );
                            }

                        /* re-set the counter */
                        pcount = 0;

                        }
                #endif
                    
                }
        
            }
    
        } /* loop over all particles */
        
    /* are there any leftovers? */
    #if VEC2
        if ( pcount > 0 ) {

            /* copy the first potential to the last entries */
            for ( k = pcount ; k < 8 ; k++ )
                potv[k] = potv[0];

            /* evaluate the potentials */
            potential_eval_vec2( potv , r2v , ev , fv );

            /* for each entry, update the forces */
            for ( k = 0 ; k < pcount ; k++ ) {
                tempv = spu_mul( dxv[k] , vec_splat( fv[k/4] , k%4 ) );
                *effi[k] = spu_sub( *effi[k] , tempv );
                *effj[k] = spu_add( *effj[k] , tempv );
                }

            }
    #else
        if ( pcount > 0 ) {

            /* copy the first potential to the last entries */
            for ( k = pcount ; k < 4 ; k++ )
                potv[k] = potv[0];

            /* evaluate the potentials */
            potential_eval_vec( potv , r2v , &ev , &fv );

            /* for each entry, update the forces */
            for ( k = 0 ; k < pcount ; k++ ) {
                tempv = spu_mul( dxv[k] , vec_splat( fv , k ) );
                *effi[k] = spu_sub( *effi[k] , tempv );
                *effj[k] = spu_add( *effj[k] , tempv );
                }

            }
    #endif
        
    }


/*////////////////////////////////////////////////////////////////////////////// */
/* int dopair_vec */
//
/* compute the interactions between a pair of cells */
/*////////////////////////////////////////////////////////////////////////////// */

void dopair ( int ni , int nj , struct part *pi , struct part *pj , float *shift ) {

    int i, j, k, count = 0, pioff, emt;
    float cutoff2 = data->cutoff * data->cutoff;
    float r2;
    #ifdef VEC2
        vector float *effi[8], *effj[8], dx, dxv[8], fv[2], ev[2], r2v[2], shiftv, tempv, pix;
    #else
        vector float *effi[4], *effj[4], dx, dxv[4], fv, ev, r2v, shiftv, tempv, pix;
    #endif
    struct part *part_i, *part_j;
    struct potential *pot, *potv[8];
    
    /* get some useful data */
    emt = data->max_type;

    /* is this a genuine pair or a cell against itself */
    if ( pj == NULL ) {
    
        /* loop over all particles */
        for ( i = 1 ; i < ni ; i++ ) {
        
            /* get the particle */
            part_i = &( pi[i] );
            pix = part_i->x;
            pioff = part_i->type * emt;
        
            /* loop over all other particles */
            for ( j = 0 ; j < i ; j++ ) {
            
                /* get the other particle */
                part_j = &( pi[j] );
                
                /* get the distance between both particles */
                dx = spu_sub( pix , part_j->x );
                r2 = _dot_product3( dx , dx );
                
                /* is this within cutoff? */
                if ( __builtin_expect( r2 > cutoff2 , 0 ) )
                    continue;
                
                /* fetch the potential, if any */
                pot = data->p[ pioff + part_j->type ];
                if ( __builtin_expect( pot == NULL , 0 ) )
                    continue;
                    
                /* add this interaction to the interaction queue */
                #ifdef VEC2
                    r2v[count/4][count%4] = r2;
                    dxv[count] = dx;
                    effi[count] = &( part_i->f );
                    effj[count] = &( part_j->f );
                    potv[count] = pot;
                    count += 1;
                    // rcount += 1;

                    /* do we have a full set to evaluate? */
                    if ( __builtin_expect( count == 8 , 0 ) ) {

                        /* evaluate the potentials */
                        potential_eval_vec2( potv , r2v , ev , fv );

                        /* for each entry, update the forces */
                        for ( k = 0 ; k < 8 ; k++ ) {
                            tempv = spu_mul( dxv[k] , vec_splat( fv[k/4] , k%4 ) );
                            *effi[k] = spu_sub( *effi[k] , tempv );
                            *effj[k] = spu_add( *effj[k] , tempv );
                            }

                        /* re-set the counter */
                        count = 0;

                        }
                #else
                    r2v[count] = r2;
                    dxv[count] = dx;
                    effi[count] = &( part_i->f );
                    effj[count] = &( part_j->f );
                    potv[count] = pot;
                    count += 1;
                    // rcount += 1;

                    /* do we have a full set to evaluate? */
                    if ( __builtin_expect( count == 4 , 0 ) ) {

                        /* evaluate the potentials */
                        potential_eval_vec( potv , r2v , &ev , &fv );

                        /* for each entry, update the forces */
                        for ( k = 0 ; k < 4 ; k++ ) {
                            tempv = spu_mul( dxv[k] , vec_splat( fv , k ) );
                            *effi[k] = spu_sub( *effi[k] , tempv );
                            *effj[k] = spu_add( *effj[k] , tempv );
                            }

                        /* re-set the counter */
                        count = 0;

                        }
                #endif
                    
                } /* loop over all other particles */
        
            } /* loop over all particles */
    
        }
        
    /* no, it's a genuine pair */
    else {
    
        /* extract shiftv */
        shiftv = spu_splats( 0.0f );
        shiftv = spu_insert( shift[0] , shiftv , 0 );
        shiftv = spu_insert( shift[1] , shiftv , 1 );
        shiftv = spu_insert( shift[2] , shiftv , 2 );
        
        /* init r2v and make the compiler happy */
        #ifdef VEC2
            r2v[0] = spu_splats( 0.0f );
            r2v[1] = spu_splats( 0.0f );
        #else
            r2v = spu_splats( 0.0f );
        #endif
    
        /* loop over all particles */
        for ( i = 0 ; i < ni ; i++ ) {
        
            /* get the particle */
            part_i = &( pi[i] );
            pix = spu_sub( part_i->x , shiftv );
            pioff = part_i->type * emt;
        
            /* loop over all other particles */
            for ( j = 0 ; j < nj ; j++ ) {
            
                /* get the other particle */
                part_j = &( pj[j] );
                
                /* get the distance between both particles */
                dx = spu_sub( pix , part_j->x );
                r2 = _dot_product3( dx , dx );
                
                /* is this within cutoff? */
                if ( __builtin_expect( r2 > cutoff2 , 1 ) )
                    continue;
                
                /* fetch the potential, if any */
                pot = data->p[ pioff + part_j->type ];
                if ( __builtin_expect( pot == NULL , 0 ) )
                    continue;
                    
                /* add this interaction to the interaction queue */
                #ifdef VEC2
                    r2v[count/4][count%4] = r2;
                    dxv[count] = dx;
                    effi[count] = &( part_i->f );
                    effj[count] = &( part_j->f );
                    potv[count] = pot;
                    count += 1;
                    rcount += 1;

                    /* do we have a full set to evaluate? */
                    if ( __builtin_expect( count == 8 , 0 ) ) {

                        /* evaluate the potentials */
                        potential_eval_vec2( potv , r2v , ev , fv );

                        /* for each entry, update the forces */
                        for ( k = 0 ; k < 8 ; k++ ) {
                            tempv = spu_mul( dxv[k] , vec_splat( fv[k/4] , k%4 ) );
                            *effi[k] = spu_sub( *effi[k] , tempv );
                            *effj[k] = spu_add( *effj[k] , tempv );
                            }

                        /* re-set the counter */
                        count = 0;

                        }
                #else
                    r2v[count] = r2;
                    dxv[count] = dx;
                    effi[count] = &( part_i->f );
                    effj[count] = &( part_j->f );
                    potv[count] = pot;
                    count += 1;
                    rcount += 1;

                    /* do we have a full set to evaluate? */
                    if ( __builtin_expect( count == 4 , 0 ) ) {

                        /* evaluate the potentials */
                        potential_eval_vec( potv , r2v , &ev , &fv );

                        /* for each entry, update the forces */
                        for ( k = 0 ; k < 4 ; k++ ) {
                            tempv = spu_mul( dxv[k] , vec_splat( fv , k ) );
                            *effi[k] = spu_sub( *effi[k] , tempv );
                            *effj[k] = spu_add( *effj[k] , tempv );
                            }

                        /* re-set the counter */
                        count = 0;

                        }
                #endif
                    
                } /* loop over all other particles */
        
            } /* loop over all particles */

        }
        
    /* are there any leftovers? */
    #ifdef VEC2
    if ( count > 0 ) {
    
        /* copy the first potential to the last entries */
        for ( k = count ; k < 8 ; k++ )
            potv[k] = potv[0];
            
        /* evaluate the potentials */
        potential_eval_vec2( potv , r2v , ev , fv );

        /* for each entry, update the forces */
        for ( k = 0 ; k < count ; k++ ) {
            tempv = spu_mul( dxv[k] , vec_splat( fv[k/4] , k%4 ) );
            *effi[k] = spu_sub( *effi[k] , tempv );
            *effj[k] = spu_add( *effj[k] , tempv );
            }
            
        }
    #else
    if ( count > 0 ) {
    
        /* copy the first potential to the last entries */
        for ( k = count ; k < 4 ; k++ )
            potv[k] = potv[0];
            
        /* evaluate the potentials */
        potential_eval_vec( potv , r2v , &ev , &fv );

        /* for each entry, update the forces */
        for ( k = 0 ; k < count ; k++ ) {
            tempv = spu_mul( dxv[k] , vec_splat( fv , k ) );
            *effi[k] = spu_sub( *effi[k] , tempv );
            *effj[k] = spu_add( *effj[k] , tempv );
            }
            
        }
    #endif
        
    }

#else
/*////////////////////////////////////////////////////////////////////////////// */
/* sortedpair */
//
/* compute the pairwise interactions for the given pair using the sorted */
/* interactions algorithm */
/*////////////////////////////////////////////////////////////////////////////// */

void sortedpair ( int ni , int nj , struct part *pi , struct part *pj , float *pshift ) {

    struct part *part_i, *part_j;
    struct potential *pot;
    float cutoff, cutoff2;
    #ifdef USE_DOUBLES
        double r2, dx[3], e, f, pjx[3];
    #else
        float r2, dx[3], e, f, pjx[3];
    #endif
    float d[2*maxparts], temp, pivot;
    float shift[3];
    int ind[2*maxparts], left[2*maxparts], count = 0, lcount = 0;
    int i, j, k, imax, qpos, lo, hi, pjoff, emt;
    struct {
        int lo, hi;
        } qstack[maxqstack];
    
    /* get the space and cutoff */
    cutoff = data->cutoff;
    cutoff2 = cutoff * cutoff;
    emt = data->max_type;
        
    /* start by filling the particle ids of both cells into ind and d */
    temp = 1.0 / sqrt( pshift[0]*pshift[0] + pshift[1]*pshift[1] + pshift[2]*pshift[2] );
    shift[0] = pshift[0]*temp; shift[1] = pshift[1]*temp; shift[2] = pshift[2]*temp;
    for ( i = 0 ; i < ni ; i++ ) {
        part_i = &( pi[i] );
        ind[count] = -i - 1;
        d[count] = part_i->x[0]*shift[0] + part_i->x[1]*shift[1] + part_i->x[2]*shift[2];
        count += 1;
        }
    for ( i = 0 ; i < nj ; i++ ) {
        part_i = &( pj[i] );
        ind[count] = i;
        d[count] = (part_i->x[0]+pshift[0])*shift[0] + (part_i->x[1]+pshift[1])*shift[1] + (part_i->x[2]+pshift[2])*shift[2] - cutoff;
        count += 1;
        }
        
    /* sort with quicksort */
    qstack[0].lo = 0; qstack[0].hi = count - 1; qpos = 0;
    while ( qpos >= 0 ) {
        lo = qstack[qpos].lo; hi = qstack[qpos].hi;
        qpos -= 1;
        if ( hi - lo < 10 ) {
            for ( i = lo ; i < hi ; i++ ) {
                imax = i;
                for ( j = i+1 ; j <= hi ; j++ )
                    if ( d[j] > d[imax] )
                        imax = j;
                if ( imax != i ) {
                    k = ind[imax]; ind[imax] = ind[i]; ind[i] = k;
                    temp = d[imax]; d[imax] = d[i]; d[i] = temp;
                    }
                }
            }
        else {
            pivot = d[ ( lo + hi ) / 2 ];
            i = lo; j = hi;
            while ( i <= j ) {
                while ( d[i] > pivot ) i++;
                while ( d[j] < pivot ) j--;
                if ( i <= j ) {
                    if ( i < j ) {
                        k = ind[i]; ind[i] = ind[j]; ind[j] = k;
                        temp = d[i]; d[i] = d[j]; d[j] = temp;
                        }
                    i += 1; j -= 1;
                    }
                }
            if ( lo < j ) {
                qpos += 1;
                qstack[qpos].lo = lo;
                qstack[qpos].hi = j;
                }
            if ( i < hi ) {
                qpos += 1;
                qstack[qpos].lo = i;
                qstack[qpos].hi = hi;
                }
            }
        }
        
    /* sort with selection sort */
    /* for ( i = 0 ; i < count-1 ; i++ ) { */
    /*     imax = i; */
    /*     for ( j = i+1 ; j < count ; j++ ) */
    /*         if ( d[j] > d[imax] ) */
    /*             imax = j; */
    /*     if ( imax != i ) { */
    /*         k = ind[imax]; ind[imax] = ind[i]; ind[i] = k; */
    /*         temp = d[imax]; d[imax] = d[i]; d[i] = temp; */
    /*         } */
    /*     } */
        
    /* loop over the sorted list of particles */
    for ( i = 0 ; i < count ; i++ ) {
    
        /* is this a particle from the left? */
        if ( ind[i] < 0 )
            left[lcount++] = -ind[i] - 1;
            
        /* it's from the right, interact with all left particles */
        else {
        
            /* get a handle on this particle */
            part_j = &( pj[ind[i]] );
            pjx[0] = part_j->x[0] + pshift[0];
            pjx[1] = part_j->x[1] + pshift[1];
            pjx[2] = part_j->x[2] + pshift[2];
            pjoff = part_j->type * emt;
        
            /* loop over the left particles */
            for ( j = lcount-1 ; j >= 0 ; j-- ) {
            
                /* get a handle on the second particle */
                part_i = &( pi[left[j]] );
            
                /* get the distance between both particles */
                for ( r2 = 0.0 , k = 0 ; k < 3 ; k++ ) {
                    dx[k] = part_i->x[k] - pjx[k];
                    r2 += dx[k] * dx[k];
                    }
                    
                /* is this within cutoff? */
                if ( r2 >= cutoff2 )
                    continue;
                
                /* fetch the potential, if any */
                pot = data->p[ pjoff + part_i->type ];
                if ( pot == NULL )
                    continue;
                    
                /* evaluate the interaction */
                /* rcount += 1; */
                #ifdef EXPLICIT_POTENTIALS
                    potential_eval_expl( pot , r2 , &e , &f );
                #else
                    potential_eval( pot , r2 , &e , &f );
                #endif
                
                /* update the forces */
                for ( k = 0 ; k < 3 ; k++ ) {
                    part_i->f[k] += -f * dx[k];
                    part_j->f[k] += f * dx[k];
                    }
                    
                }
        
            }
    
        } /* loop over all particles */
        
    }


/*////////////////////////////////////////////////////////////////////////////// */
/* int dopair */
//
/* compute the interactions between a pair of cells */
/*////////////////////////////////////////////////////////////////////////////// */

void dopair ( int ni , int nj , struct part *pi , struct part *pj , float *shift ) {

    int i, j, k, pioff, emt;
    float cutoff2 = data->cutoff * data->cutoff;
    #ifdef USE_DOUBLES
        double dx[3], r2, e, f, pix[3];
    #else
        float dx[3], r2, e, f, pix[3];
    #endif
    struct part *part_i, *part_j;
    struct potential *pot;
    
    /* get the max type */
    emt = data->max_type;

    /* is this a genuine pair or a cell against itself */
    if ( pj == NULL ) {
    
        /* loop over all particles */
        for ( i = 1 ; i < ni ; i++ ) {
        
            /* get the particle */
            part_i = &(pi[i]);
            pix[0] = part_i->x[0];
            pix[1] = part_i->x[1];
            pix[2] = part_i->x[2];
            pioff = part_i->type * emt;
        
            /* loop over all other particles */
            for ( j = 0 ; j < i ; j++ ) {
            
                /* get the other particle */
                part_j = &(pi[j]);
                
                /* get the distance between both particles */
                for ( r2 = 0.0 , k = 0 ; k < 3 ; k++ ) {
                    dx[k] = pix[k] - part_j->x[k];
                    r2 += dx[k] * dx[k];
                    }
                    
                /* is this within cutoff? */
                if ( r2 >= cutoff2 )
                    continue;
                
                /* fetch the potential, if any */
                pot = data->p[ pioff + part_j->type ];
                if ( pot == NULL )
                    continue;
                    
                /* evaluate the interaction */
                /* rcount += 1; */
                #ifdef EXPLICIT_POTENTIALS
                    potential_eval_expl( pot , r2 , &e , &f );
                #else
                    potential_eval( pot , r2 , &e , &f );
                #endif
                
                /* update the forces */
                for ( k = 0 ; k < 3 ; k++ ) {
                    part_i->f[k] += -f * dx[k];
                    part_j->f[k] += f * dx[k];
                    }
                    
                } /* loop over all other particles */
        
            } /* loop over all particles */
    
        }
        
    /* no, it's a genuine pair */
    else {
    
        /* loop over all particles */
        for ( i = 0 ; i < ni ; i++ ) {
        
            /* get the particle */
            part_i = &(pi[i]);
            pix[0] = part_i->x[0] - shift[0];
            pix[1] = part_i->x[1] - shift[1];
            pix[2] = part_i->x[2] - shift[2];
            pioff = part_i->type * emt;
        
            /* loop over all other particles */
            for ( j = 0 ; j < nj ; j++ ) {
            
                /* get the other particle */
                part_j = &(pj[j]);

                /* get the distance between both particles */
                for ( r2 = 0.0 , k = 0 ; k < 3 ; k++ ) {
                    dx[k] = pix[k] - part_j->x[k];
                    r2 += dx[k] * dx[k];
                    }
                    
                /* is this within cutoff? */
                if ( r2 >= cutoff2 )
                    continue;
                
                /* fetch the potential, if any */
                pot = data->p[ pioff + part_j->type ];
                if ( pot == NULL )
                    continue;
                    
                /* evaluate the interaction */
                /* rcount += 1; */
                #ifdef EXPLICIT_POTENTIALS
                    potential_eval_expl( pot , r2 , &e , &f );
                #else
                    potential_eval( pot , r2 , &e , &f );
                #endif
                
                /* update the forces */
                for ( k = 0 ; k < 3 ; k++ ) {
                    part_i->f[k] += -f * dx[k];
                    part_j->f[k] += f * dx[k];
                    }
                
                } /* loop over all other particles */
        
            } /* loop over all particles */

        }
        
    }
#endif

/*////////////////////////////////////////////////////////////////////////////// */
/* void my_get */
//
/* wrapper for mfc_get that checks for sizes above 16k and breaks things */
/* down to multiple calls. */
/*////////////////////////////////////////////////////////////////////////////// */

inline void my_get ( void *to , unsigned long long from , unsigned int size , unsigned int tid ) {

    /* cut off chunks while the data is too large... */
    while ( __builtin_expect( size > 0x4000 , 0 ) ) {
        mfc_get( to , from , 0x4000 , tid , 0 , 0 );
        from += 0x4000;
        to += 0x4000;
        size -= 0x4000;
        }
        
    /* send the leftovers */
    mfc_get( to , from , size , tid , 0 , 0 );

    }
    

/*////////////////////////////////////////////////////////////////////////////// */
/* void my_put */
//
/* wrapper for mfc_put that checks for sizes above 16k and breaks things */
/* down to multiple calls. */
/*////////////////////////////////////////////////////////////////////////////// */

inline void my_put ( void *from , unsigned long long to , unsigned int size , unsigned int tid ) {

    /* cut off chunks while the data is too large... */
    while ( __builtin_expect( size > 0x4000 , 0 ) ) {
        mfc_put( from , to , 0x4000 , tid , 0 , 0 );
        from += 0x4000;
        to += 0x4000;
        size -= 0x4000;
        }
        
    /* send the leftovers */
    mfc_put( from , to , size , tid , 0 , 0 );

    }
    

/*////////////////////////////////////////////////////////////////////////////// */
/* int main */
//
/* this is the main entry point for the spu routine */
/*////////////////////////////////////////////////////////////////////////////// */

int main ( unsigned long long id , unsigned long long argp , unsigned long long envp ) {

    unsigned int tid[bitesize];
    int size_data = envp;
    int i, j, k, bid = 0, nid = 1, nnid = 2;
    int ni[bitesize], nj[bitesize], status[bitesize], wbi[bitesize], wbj[bitesize];
    unsigned long long ai[bitesize], aj[bitesize];
    struct part *pi[bitesize], *pj[bitesize], *pbuff[2*bitesize];
    float shift[3*bitesize];
    unsigned int buff, dummy;
    int pbuff_first = 0, pbuff_last = 0;
    struct cell *cells;

    /* say hello */
    printf("runner_spu: spu 0x%llx says hi.\n",id); fflush(stdout);
    #ifdef VECTORIZE
        printf("runner_spu: VECTORIZE is on.\n");
    #endif
    #ifdef SORTED_INTERACTIONS
        printf("runner_spu: SORTED_INTERACTIONS is on.\n");
    #endif
    
    
    /* allocate some local space for the engine data */
    if ( ( data = (struct data *)malloc_align( size_data , 7 ) ) == NULL )
        return -__LINE__;
        
    /* reserve the DMA tags */
    for ( k = 0 ; k < bitesize ; k++ )
        if ( ( tid[k] = mfc_tag_reserve() ) == MFC_TAG_INVALID )
            return -__LINE__;
    
    /* make a call to get the data */
    my_get( data , argp , size_data , tid[0] );
    mfc_write_tag_mask( 1 << tid[0] );
    mfc_read_tag_status_all();
    
    /* say something */
    /* printf("runner_spu: loaded %i bytes of data from main store.\n", size_data); fflush(stdout); */
    
    /* hook-up the data structures */
    for ( i = 0 ; i < data->max_type * data->max_type ; i++ )
        if ( data->p[i] != NULL )
            data->p[i] = (struct potential *)( (unsigned int)data + (unsigned int)data->p[i] );
        else
            data->p[i] = NULL;
            
    /* allocate space for the cells list */
    if ( ( cells = (struct cell *)malloc_align( mfc_ceil128( sizeof(struct cell) * data->nr_cells ) , 7 ) ) == NULL )
        return -__LINE__;
            
            
    /* init some stuff */
    for ( k = 0 ; k < bitesize ; k++ ) {
        ai[k] = 0; aj[k] = 0;
        pi[k] = NULL; pj[k] = NULL;
        status[k] = 0;
        }
    for ( k = 0 ; k < 2*bitesize ; k++ )
        if ( ( pbuff[k] = (struct part *)malloc_align( mfc_ceil128( sizeof(struct part) * maxparts ) , 7 ) ) == NULL )
            return -__LINE__;
        
    /* tell the PPU that we are ready. */
    spu_write_out_intr_mbox( 1 );
    
    
            
    /* this is the main loop which picks up a pair and computes its */
    /* interactions. */
    while ( 1 ) {
    
        /* First step: feed third pipe. */
        
        /* check if there is new data waiting */
        if ( ( spu_stat_in_mbox() > 0 ) ||
            ( status[bid] == 0 && status[nid] == 0 && status[nnid] == 0 ) ) {
            
            /* get the data for the next pair */
            /* printf("runner_spu: runner 0x%llx waiting on mailbox input...\n",id); fflush(stdout); */
            buff = spu_read_in_mbox();
            dummy = spu_read_in_mbox();
            /* printf("runner_spu: runner 0x%llx got mbox input.\n",id); fflush(stdout); */

            /* did we get a flush message? */
            if ( buff == 0 ) {
            
                /* mark this... */
                /* printf("runner_run: spu 0x%llx got a flush message...\n",id); fflush(stdout); */
                status[nnid] = -1;
                
                }
                
            /* did we get a reload message? */
            else if ( buff == 0xFFFFFFFF ) {
            
                /* printf("runner_run: spu 0x%llx got a reload message...\n",id); fflush(stdout); */
                
                /* emit mfc */
                my_get( cells , data->cell_addr , mfc_ceil128( sizeof(struct cell) * data->nr_cells ) , tid[nnid] );
                
                /* wait for the data */
                mfc_write_tag_mask( 1 << tid[nnid] );
                mfc_read_tag_status_all();
                
                /* carry on as if nothing happened... */
                /* printf("runner_run: spu 0x%llx executed reload.\n",id); fflush(stdout); */
                continue;
            
                }
                
            /* we got a genuine cell pair  */
            else {

                /* unpack the data */
                i = buff >> 20;
                j = (buff >> 8) & 4095;
                ni[nnid] = cells[i].ni;
                nj[nnid] = cells[j].ni;
                ai[nnid] = cells[i].ai;
                aj[nnid] = cells[j].ai;
                if ( ( (buff >> 6) & 3 ) == 0 )
                    shift[nnid*3+0] = 0.0;
                else if ( ( (buff >> 6) & 3 ) == 1 )
                    shift[nnid*3+0] = data->h[0];
                else
                    shift[nnid*3+0] = -data->h[0];
                if ( ( (buff >> 4) & 3 ) == 0 )
                    shift[nnid*3+1] = 0.0;
                else if ( ( (buff >> 4) & 3 ) == 1 )
                    shift[nnid*3+1] = data->h[1];
                else
                    shift[nnid*3+1] = -data->h[1];
                if ( ( (buff >> 2) & 3 ) == 0 )
                    shift[nnid*3+2] = 0.0;
                else if ( ( (buff >> 2) & 3 ) == 1 )
                    shift[nnid*3+2] = data->h[2];
                else
                    shift[nnid*3+2] = -data->h[2];
                /* printf("runner_spu: got pair 0x%llx (n=%i), 0x%llx (n=%i) with shift=[%e,%e,%e].\n",
                    ai[nnid],ni[nnid],aj[nnid],nj[nnid],shift[nnid*3+0],shift[nnid*3+1],shift[nnid*3+2]); fflush(stdout); */

                /* check if we can recycle any of these buffers */
                pi[nnid] = NULL; pj[nnid] = NULL;
                wbi[nnid] = 1; wbj[nnid] = ( ai[nnid] != aj[nnid] );
                if ( status[nid] > 0 && ai[nid] == ai[nnid] ) {
                    pi[nnid] = pi[nid]; wbi[nid] = 0;
                    }
                else if ( status[nid] > 0 && aj[nid] == ai[nnid] ) {
                    pi[nnid] = pj[nid]; wbj[nid] = 0;
                    }
                else if ( status[bid] > 0 && ai[bid] == ai[nnid] ) {
                    pi[nnid] = pi[bid]; wbi[bid] = 0;
                    }
                else if ( status[bid] > 0 && aj[bid] == ai[nnid] ) {
                    pi[nnid] = pj[bid]; wbj[bid] = 0;
                    }
                if ( wbj[nnid] ) {
                    if ( status[nid] > 0 && ai[nid] == aj[nnid] ) {
                        pj[nnid] = pi[nid]; wbi[nid] = 0;
                        }
                    else if ( status[nid] > 0 && aj[nid] == aj[nnid] ) {
                        pj[nnid] = pj[nid]; wbj[nid] = 0;
                        }
                    else if ( status[bid] > 0 && ai[bid] == aj[nnid] ) {
                        pj[nnid] = pi[bid]; wbi[bid] = 0;
                        }
                    else if ( status[bid] > 0 && aj[bid] == aj[nnid] ) {
                        pj[nnid] = pj[bid]; wbj[bid] = 0;
                        }
                    }

                /* emit the load for the particles of this pair if needed */
                if ( pi[nnid] == NULL ) {
                    pi[nnid] = pbuff[ pbuff_first ]; pbuff_first = ( pbuff_first + 1 ) % ( 2 * bitesize );
                    my_get( pi[nnid] , ai[nnid] , mfc_ceil128( sizeof(struct part) * ni[nnid] ) , tid[nnid] );
                    }
                if ( wbj[nnid] && pj[nnid] == NULL ) {
                    pj[nnid] = pbuff[ pbuff_first ]; pbuff_first = ( pbuff_first + 1 ) % ( 2 * bitesize );
                    my_get( pj[nnid] , aj[nnid] , mfc_ceil128( sizeof(struct part) * nj[nnid] ) , tid[nnid] );
                    }

                /* set the status */
                status[nnid] = 1;
                
                }
                
            }
                        
            
        /* Second step: wait for get, compute, and send put on second pipe */
        
        if ( status[nid] > 0 ) {
        
            /* wait for the current data */
            mfc_write_tag_mask( 1 << tid[nid] );
            mfc_read_tag_status_all();
            
            /* compute the interactions in this pair */
            if ( ai[nid] == aj[nid] )
                dopair( ni[nid] , 0 , pi[nid] , NULL , NULL );
            else
                #ifdef SORTED_INTERACTIONS
                    sortedpair( ni[nid] , nj[nid] , pi[nid] , pj[nid] , &( shift[nid*3] ) );
                #else
                    dopair( ni[nid] , nj[nid] , pi[nid] , pj[nid] , &( shift[nid*3] ) );
                #endif
                
            /* write the data back to the main memory and free the buffers, */
            /* if needed. */
            if ( wbi[nid] )
                my_put( pi[nid] , ai[nid] , mfc_ceil128( ni[nid] * sizeof(struct part) ) , tid[nid] );
            if ( wbj[nid] )
                my_put( pj[nid] , aj[nid] , mfc_ceil128( nj[nid] * sizeof(struct part) ) , tid[nid] );
                
            }
            
            
        /* Third step: wait for put and signal on first pipe */
        
        if ( status[bid] > 0 ) {
        
            /* wait for any data on this channel to have been written... */
            mfc_write_tag_mask( 1 << tid[bid] );
            mfc_read_tag_status_all();
            
            /* return the buffers to the queue */
            if ( wbi[bid] ) {
                pbuff[ pbuff_last ] = pi[bid];
                pbuff_last = ( pbuff_last + 1 ) % ( 2 * bitesize );
                }
            if ( wbj[bid] ) {
                pbuff[ pbuff_last ] = pj[bid];
                pbuff_last = ( pbuff_last + 1 ) % ( 2 * bitesize );
                }
            
            /* signal that we are done with this pair (if there was anything there) */
            /* printf("runner_spu: runner 0x%llx waiting on mailbox output...\n",id); fflush(stdout); */
            /* spu_write_out_intr_mbox( rcount ); */
            /* printf("runner_spu: runner 0x%llx got mailbox output...\n",id); fflush(stdout); */
            
            /* done with this box */
            status[bid] = 0;

            }
            
        /* or did we get a flush? */
        else if ( status[bid] < 0 ) {
        
            /* signal to the PPU that we are through... */
            /* printf("runner_run: spu 0x%llx executing flush...\n",id); fflush(stdout); */
            spu_write_out_intr_mbox( rcount );
            /* printf("runner_run: spu 0x%llx executed flush.\n",id); fflush(stdout); */
            
            /* done with this box */
            status[bid] = 0;
        
            }
        

        /* move on... */
        bid = nid;
        nid = nnid;
        nnid = ( nid + 1 ) % bitesize;

        }
    
    }
