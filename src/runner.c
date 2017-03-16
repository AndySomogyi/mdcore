/*******************************************************************************
 * This file is part of mdcore.
 * Coypright (c) 2010 Pedro Gonnet (pedro.gonnet@durham.ac.uk)
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

/* Include configuration header */
#include "../config.h"

/* Include some standard header files */
#include <stdlib.h>
#include <stdio.h>
#include <pthread.h>
#include <math.h>
#include <float.h>
#include <string.h>
#include <limits.h>

/* Include some conditional headers. */
#include "../config.h"
#ifdef CELL
    #include <libspe2.h>
    #include <libmisc.h>
    #define ceil128(v) (((v) + 127) & ~127)
#endif
#ifdef __SSE__
    #include <xmmintrin.h>
#endif
#ifdef HAVE_SETAFFINITY
    #include <sched.h>
#endif
#ifdef WITH_MPI
    #include <mpi.h>
#endif

/* Include local headers */
#include "cycle.h"
#include "errs.h"
#include "fptype.h"
#include "lock.h"
#include "part.h"
#include "queue.h"
#include "cell.h"
#include "task.h"
#include "space.h"
#include "potential.h"
#include "engine.h"
#include "runner.h"



#ifdef CELL
    /* the SPU executeable */
    extern spe_program_handle_t runner_spu;
#endif


/* Global variables. */
/** The ID of the last error. */
int runner_err = runner_err_ok;
unsigned int runner_rcount = 0;

/** Timers. */
ticks runner_timers[runner_timer_count];

/* the error macro. */
#define error(id)				( runner_err = errs_register( id , runner_err_msg[-(id)] , __LINE__ , __FUNCTION__ , __FILE__ ) )

/* list of error messages. */
char *runner_err_msg[12] = {
	"Nothing bad happened.",
    "An unexpected NULL pointer was encountered.",
    "A call to malloc failed, probably due to insufficient memory.",
    "An error occured when calling a space function.",
    "A call to a pthread routine failed.",
    "An error occured when calling an engine function.",
    "An error occured when calling an SPE function.",
    "An error occured with the memory flow controler.",
    "The requested functionality is not available." ,
    "An error occured when calling an fifo function." ,
    "Error filling Verlet list: too many neighbours." , 
    "Unknown task type." , 
	};
    
    

/**
 * @brief Sort the particles in descending order using QuickSort.
 *
 * @param parts The particle IDs and distances in compact form
 * @param N The number of particles.
 *
 * The particle data is assumed to contain the distance in the lower
 * 16 bits and the particle ID in the upper 16 bits.
 */
 
void runner_sort_descending ( unsigned int *parts , int N ) {

    struct {
        short int lo, hi;
        } qstack[10];
    int qpos, i, j, lo, hi, pivot, imax;
    unsigned int temp;
        
    /* Sort parts in cell_i in decreasing order with quicksort */
    qstack[0].lo = 0; qstack[0].hi = N - 1; qpos = 0;
    while ( qpos >= 0 ) {
        lo = qstack[qpos].lo; hi = qstack[qpos].hi;
        qpos -= 1;
        if ( hi - lo < 15 ) {
            for ( i = lo ; i < hi ; i++ ) {
                imax = i;
                for ( j = i+1 ; j <= hi ; j++ )
                    if ( (parts[j] & 0xffff) > (parts[imax] & 0xffff) )
                        imax = j;
                if ( imax != i ) {
                    temp = parts[imax]; parts[imax] = parts[i]; parts[i] = temp;
                    }
                }
            }
        else {
            pivot = parts[ ( lo + hi ) / 2 ] & 0xffff;
            i = lo; j = hi;
            while ( i <= j ) {
                while ( (parts[i] & 0xffff) > pivot ) i++;
                while ( (parts[j] & 0xffff) < pivot ) j--;
                if ( i <= j ) {
                    if ( i < j ) {
                        temp = parts[i]; parts[i] = parts[j]; parts[j] = temp;
                        }
                    i += 1; j -= 1;
                    }
                }
            if ( j > ( lo + hi ) / 2 ) {
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
            else {
                if ( i < hi ) {
                    qpos += 1;
                    qstack[qpos].lo = i;
                    qstack[qpos].hi = hi;
                    }
                if ( lo < j ) {
                    qpos += 1;
                    qstack[qpos].lo = lo;
                    qstack[qpos].hi = j;
                    }
                }
            }
        }
                
    }
    

/**
 * @brief Sort the particles in ascending order using QuickSort.
 *
 * @param parts The particle IDs and distances in compact form
 * @param N The number of particles.
 *
 * The particle data is assumed to contain the distance in the lower
 * 16 bits and the particle ID in the upper 16 bits.
 */
 
void runner_sort_ascending ( unsigned int *parts , int N ) {

    struct {
        short int lo, hi;
        } qstack[10];
    int qpos, i, j, lo, hi, pivot, imax;
    unsigned int temp;
        
    /* Sort parts in cell_i in decreasing order with quicksort */
    qstack[0].lo = 0; qstack[0].hi = N - 1; qpos = 0;
    while ( qpos >= 0 ) {
        lo = qstack[qpos].lo; hi = qstack[qpos].hi;
        qpos -= 1;
        if ( hi - lo < 15 ) {
            for ( i = lo ; i < hi ; i++ ) {
                imax = i;
                for ( j = i+1 ; j <= hi ; j++ )
                    if ( (parts[j] & 0xffff) < (parts[imax] & 0xffff) )
                        imax = j;
                if ( imax != i ) {
                    temp = parts[imax]; parts[imax] = parts[i]; parts[i] = temp;
                    }
                }
            }
        else {
            pivot = parts[ ( lo + hi ) / 2 ] & 0xffff;
            i = lo; j = hi;
            while ( i <= j ) {
                while ( (parts[i] & 0xffff) < pivot ) i++;
                while ( (parts[j] & 0xffff) > pivot ) j--;
                if ( i <= j ) {
                    if ( i < j ) {
                        temp = parts[i]; parts[i] = parts[j]; parts[j] = temp;
                        }
                    i += 1; j -= 1;
                    }
                }
            if ( j > ( lo + hi ) / 2 ) {
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
            else {
                if ( i < hi ) {
                    qpos += 1;
                    qstack[qpos].lo = i;
                    qstack[qpos].hi = hi;
                    }
                if ( lo < j ) {
                    qpos += 1;
                    qstack[qpos].lo = lo;
                    qstack[qpos].hi = j;
                    }
                }
            }
        }
                
    }
    

/**
 * @brief The #runner's main routine (for the Cell/BE SPU).
 *
 * @param r Pointer to the #runner to run.
 *
 * @return #runner_err_ok or <0 on error (see #runner_err).
 *
 * This is the main routine for the #runner. When called, it enters
 * an infinite loop in which it waits at the #engine @c r->e barrier
 * and, once having paSSEd, calls #space_getpair until there are no pairs
 * available.
 *
 * Note that this routine is only compiled if @c CELL has been defined.
 */

int runner_run_cell ( struct runner *r ) {

#ifdef CELL
    int err = 0;
    struct cellpair *p[runner_qlen];
    unsigned int buff[2];
    int i, k, count = 0;
    struct space *s;

    /* check the inputs */
    if ( r == NULL )
        return error(runner_err_null);
        
    /* init some local pointers. */
    s = &(r->e->s);
        
    /* give a hoot */
    printf("runner_run: runner %i is up and running (SPU)...\n",r->id); fflush(stdout);
    
    /* init the cellpair pointers */
    for ( k = 0 ; k < runner_qlen ; k++ )
        p[k] = NULL;
        
    /* main loop, in which the runner should stay forever... */
    while ( 1 ) {
    
        /* wait at the engine barrier */
        /* printf("runner_run: runner %i waiting at barrier...\n",r->id); */
        if ( engine_barrier(r->e) < 0)
            return error(runner_err_engine);
            
        /* write the current cell data */
        for ( i = 0 ; i < s->nr_cells ; i++ ) {
            r->celldata[i].ni = s->cells[i].count;
            r->celldata[i].ai = (unsigned long long)s->cells[i].parts;
            }

        /* emit a reload message */
        buff[0] = 0xFFFFFFFF;
        /* printf("runner_run: runner %i sending reload message...\n",r->id); */
        if ( spe_in_mbox_write( r->spe , buff , 2 , SPE_MBOX_ALL_BLOCKING ) != 2 )
            return runner_err_spe;


        /* while there are pairs... */
        while ( s->next_pair < s->nr_pairs || count > 0 ) {

            /* if we have no p[0], try to get some... */
            if ( p[0] == NULL && s->next_pair < s->nr_pairs ) {
                p[0] = space_getpair( &(r->e->s) , r->id , runner_bitesize , NULL , &err , count == 0 );
                if ( err < 0 )
                    return runner_err_space;
                }

            /* if we got a pair, send it to the SPU... */
            if ( p[0] != NULL ) {

                /* we've got an active slot! */
                count += 1;

                /* pack this pair's data */
                buff[0] = ( p[0]->i << 20 ) + ( p[0]->j << 8 ) + 1;
                if ( p[0]->shift[0] > 0 )
                    buff[0] += 1 << 6;
                else if ( p[0]->shift[0] < 0 )
                    buff[0] += 2 << 6;
                if ( p[0]->shift[1] > 0 )
                    buff[0] += 1 << 4;
                else if ( p[0]->shift[1] < 0 )
                    buff[0] += 2 << 4;
                if ( p[0]->shift[2] > 0 )
                    buff[0] += 1 << 2;
                else if ( p[0]->shift[2] < 0 )
                    buff[0] += 2 << 2;

                /* wait for the buffer to be free... */
                /* while ( !spe_in_mbox_status( r->spe ) ) */
                /*     sched_yield(); */

                /* write the data to the mailbox */
                /* printf("runner_run: sending pair 0x%llx (n=%i), 0x%llx (n=%i) with shift=[%e,%e,%e].\n",
                    (unsigned long long)s->cells[p[0]->i].parts,s->cells[p[0]->i].count,(unsigned long long)s->cells[p[0]->j].parts,s->cells[p[0]->j].count,
                    p[0]->shift[0], p[0]->shift[1], p[0]->shift[2]); fflush(stdout); */
                /* printf("runner_run: runner %i sending pair to SPU...\n",r->id); fflush(stdout); */
                if ( spe_in_mbox_write( r->spe , buff , 2 , SPE_MBOX_ALL_BLOCKING ) != 2 )
                    return runner_err_spe;
                /* printf("runner_run: runner %i sent pair to SPU.\n",r->id); fflush(stdout); */


                /* wait for the last pair to have been proceSSEd */
                if ( p[runner_qlen-1] != NULL ) {

                    /* read a word from the spe */
                    /* printf("runner_run: runner %i waiting for SPU response...\n",r->id); fflush(stdout); */
                    /* if ( spe_out_intr_mbox_read( r->spe , &buff , 1 , SPE_MBOX_ALL_BLOCKING ) < 1 )
                        return runner_err_spe; */
                    /* printf("runner_run: runner %i got SPU response.\n",r->id); fflush(stdout); */

                    /* release the last pair */
                    if ( space_releasepair( s , p[runner_qlen-1]->i , p[runner_qlen-1]->j ) < 0 )
                        return runner_err_space;

                    /* we've got one less... */
                    count -= 1;

                    }

                /* move on in the chain */
                for ( k = runner_qlen-1 ; k > 0 ; k-- )
                    p[k] = p[k-1];
                if ( p[0] != NULL )
                    p[0] = p[0]->next;

                /* take a breather... */
                /* sched_yield(); */

                }

            /* is there a non-empy slot, send a flush */
            else if ( count > 0 ) {

                /* send a flush message... */
                buff[0] = 0;
                if ( spe_in_mbox_write( r->spe , buff , 2 , SPE_MBOX_ALL_BLOCKING ) != 2 )
                    return runner_err_spe;

                /* wait for the reply... */
                if ( spe_out_intr_mbox_read( r->spe , buff , 1 , SPE_MBOX_ALL_BLOCKING ) < 1 )
                    return runner_err_spe;
                /* printf("runner_run: got rcount=%u.\n",buff[0]); */

                /* release the pairs still in the queue */
                for ( k = 1 ; k < runner_qlen ; k++ )
                    if ( p[k] != NULL ) {
                        if ( space_releasepair( &(r->e->s) , p[k]->i , p[k]->j ) < 0 )
                            return runner_err_space;
                        p[k] = NULL;
                        count -= 1;
                        }

                }

            }
                
        /* did things go wrong? */
        /* printf("runner_run: runner %i done pairs.\n",r->id); fflush(stdout); */
        if ( err < 0 )
            return error(runner_err_space);
    
        }

    /* end well... */
    return runner_err_ok;

#else

    /* This functionality is not available */
    return runner_err_unavail;
    
#endif

    }


/**
 * @brief The #runner's main routine (for the Cell/BE SPU, using
 *      tuples).
 *
 * @param r Pointer to the #runner to run.
 *
 * @return #runner_err_ok or <0 on error (see #runner_err).
 *
 * This is the main routine for the #runner. When called, it enters
 * an infinite loop in which it waits at the #engine @c r->e barrier
 * and, once having paSSEd, calls #space_gettuple until there are no
 * tuples available.
 *
 * Note that this routine is only compiled if @c CELL has been defined.
 */

int runner_run_cell_tuples ( struct runner *r ) {

#ifdef CELL
    int err = 0;
    unsigned int buff[2];
    int i, j, k, count = 0, res;
    struct space *s;
    struct celltuple *t;
    int cid[runner_qlen][2];
    FPTYPE shift[3];
    

    /* check the inputs */
    if ( r == NULL )
        return error(runner_err_null);
        
    /* init some local pointers. */
    s = &(r->e->s);
        
    /* give a hoot */
    printf("runner_run: runner %i is up and running (SPU)...\n",r->id); fflush(stdout);
    
    /* init the cellpair pointers */
    for ( k = 0 ; k < runner_qlen ; k++ )
        cid[k][0] = -1;
        
    /* main loop, in which the runner should stay forever... */
    while ( 1 ) {
    
        /* wait at the engine barrier */
        /* printf("runner_run: runner %i waiting at barrier...\n",r->id); */
        if ( engine_barrier(r->e) < 0)
            return error(runner_err_engine);
            
        /* write the current cell data */
        for ( i = 0 ; i < s->nr_cells ; i++ ) {
            r->celldata[i].ni = s->cells[i].count;
            r->celldata[i].ai = (unsigned long long)s->cells[i].parts;
            }

        /* emit a reload message */
        buff[0] = 0xFFFFFFFF;
        /* printf("runner_run: runner %i sending reload message...\n",r->id); */
        if ( spe_in_mbox_write( r->spe , buff , 2 , SPE_MBOX_ALL_BLOCKING ) != 2 )
            return runner_err_spe;


        /* Loop over tuples... */
        while ( 1 ) {
        
            /* Try to get a tuple. */
            if ( ( res = space_gettuple( s , &t , count == 0 ) ) < 0 )
                return r->err = runner_err_space;
                
            /* Did we get a tuple back? */
            if ( res > 0 )
                
                /* Loop over the cell pairs in this tuple. */
                for ( i = 0 ; i < t->n ; i++ )
                    for ( j = i ; j < t->n ; j++ ) {

                        /* Is this pair active? */
                        if ( t->pairid[ space_pairind(i,j) ] < 0 )
                            continue;

                        /* Get the cell ids. */
                        cid[0][0] = t->cellid[i];
                        cid[0][1] = t->cellid[j];

                        /* we've got an active slot! */
                        count += 1;

                        /* Compute the shift between ci and cj. */
                        for ( k = 0 ; k < 3 ; k++ ) {
                            shift[k] = s->cells[cid[0][1]].origin[k] - s->cells[cid[0][0]].origin[k];
                            if ( shift[k] * 2 > s->dim[k] )
                                shift[k] -= s->dim[k];
                            else if ( shift[k] * 2 < -s->dim[k] )
                                shift[k] += s->dim[k];
                            }

                        /* pack this pair's data */
                        buff[0] = ( cid[0][0] << 20 ) + ( cid[0][1] << 8 ) + 1;
                        if ( shift[0] > 0 )
                            buff[0] += 1 << 6;
                        else if ( shift[0] < 0 )
                            buff[0] += 2 << 6;
                        if ( shift[1] > 0 )
                            buff[0] += 1 << 4;
                        else if ( shift[1] < 0 )
                            buff[0] += 2 << 4;
                        if ( shift[2] > 0 )
                            buff[0] += 1 << 2;
                        else if ( shift[2] < 0 )
                            buff[0] += 2 << 2;

                        /* write the data to the mailbox */
                        /* printf("runner_run: sending pair 0x%llx (n=%i), 0x%llx (n=%i) with shift=[%e,%e,%e].\n",
                            (unsigned long long)s->cells[p[0]->i].parts,s->cells[p[0]->i].count,(unsigned long long)s->cells[p[0]->j].parts,s->cells[p[0]->j].count,
                            p[0]->shift[0], p[0]->shift[1], p[0]->shift[2]); fflush(stdout); */
                        /* printf("runner_run: runner %i sending pair to SPU...\n",r->id); fflush(stdout); */
                        if ( spe_in_mbox_write( r->spe , buff , 2 , SPE_MBOX_ALL_BLOCKING ) != 2 )
                            return runner_err_spe;
                        /* printf("runner_run: runner %i sent pair to SPU.\n",r->id); fflush(stdout); */


                        /* wait for the last pair to have been processed */
                        if ( cid[runner_qlen-1][0] >= 0 ) {

                            /* read a word from the spe */
                            /* printf("runner_run: runner %i waiting for SPU response...\n",r->id); fflush(stdout); */
                            /* if ( spe_out_intr_mbox_read( r->spe , &buff , 1 , SPE_MBOX_ALL_BLOCKING ) < 1 )
                                return runner_err_spe; */
                            /* printf("runner_run: runner %i got SPU response.\n",r->id); fflush(stdout); */

                            /* release the last pair */
                            if ( space_releasepair( s , cid[runner_qlen-1][0] , cid[runner_qlen-1][1] ) < 0 )
                                return runner_err_space;

                            /* we've got one less... */
                            count -= 1;

                            }

                        /* move on in the chain */
                        for ( k = runner_qlen-1 ; k > 0 ; k-- ) {
                            cid[k][0] = cid[k-1][0];
                            cid[k][1] = cid[k-1][1];
                            }
                        cid[0][0] = -1;

                        }
            
            /* Did we get a stall? */
            else if ( s->next_tuple < s->nr_tuples ) {
            
                /* wait for the last pair to have been processed */
                if ( cid[runner_qlen-1][0] >= 0 ) {

                    /* read a word from the spe */
                    /* printf("runner_run: runner %i waiting for SPU response...\n",r->id); fflush(stdout); */
                    /* if ( spe_out_intr_mbox_read( r->spe , &buff , 1 , SPE_MBOX_ALL_BLOCKING ) < 1 )
                        return runner_err_spe; */
                    /* printf("runner_run: runner %i got SPU response.\n",r->id); fflush(stdout); */

                    /* release the last pair */
                    if ( space_releasepair( s , cid[runner_qlen-1][0] , cid[runner_qlen-1][1] ) < 0 )
                        return runner_err_space;

                    /* we've got one less... */
                    count -= 1;

                    }

                /* move on in the chain */
                for ( k = runner_qlen-1 ; k > 0 ; k-- ) {
                    cid[k][0] = cid[k-1][0];
                    cid[k][1] = cid[k-1][1];
                    }
                cid[0][0] = -1;

                }
                
            /* Otherwise, we're done. */
            else
                break;
                    
            }

        /* If there is a non-empy slot, send a flush */
        if ( count > 0 ) {

            /* send a flush message... */
            buff[0] = 0;
            if ( spe_in_mbox_write( r->spe , buff , 2 , SPE_MBOX_ALL_BLOCKING ) != 2 )
                return runner_err_spe;

            /* wait for the reply... */
            if ( spe_out_intr_mbox_read( r->spe , buff , 1 , SPE_MBOX_ALL_BLOCKING ) < 1 )
                return runner_err_spe;
            // printf("runner_run: got rcount=%u.\n",buff[0]);

            /* release the pairs still in the queue */
            for ( k = 1 ; k < runner_qlen ; k++ )
                if ( cid[k][0] >= 0 ) {
                    if ( space_releasepair( &(r->e->s) , cid[k][0] , cid[k][1] ) < 0 )
                        return runner_err_space;
                    cid[k][0] = -1;
                    count -= 1;
                    }

            }

        /* did things go wrong? */
        /* printf("runner_run: runner %i done pairs.\n",r->id); fflush(stdout); */
        if ( err < 0 )
            return error(runner_err_space);
    
        }

    /* end well... */
    return runner_err_ok;

#else

    /* This functionality is not available */
    return runner_err_unavail;
    
#endif

    }


int runner_run ( struct runner *r ) {

    struct engine *e = r->e;
    struct space *s = &e->s;
    int k, err = 0, acc = 0, naq, qid, myqid = e->nr_queues * r->id / e->nr_runners;
    struct task *t = NULL;
    struct queue *myq = &e->queues[ myqid ], *queues[ e->nr_queues ];
    unsigned int myseed = rand() + r->id;
    int count;

    /* give a hoot */
    printf( "runner_run: runner %i is up and running on queue %i (tasks)...\n" , r->id , myqid ); fflush(stdout);
    
    /* main loop, in which the runner should stay forever... */
    while ( 1 ) {
    
        /* wait at the engine barrier */
        /* printf("runner_run: runner %i waiting at barrier...\n",r->id); */
        if ( engine_barrier(e) < 0)
            return error(runner_err_engine);
            
        /* Init the reaction counter. */
        // runner_rcount = 0;
        
        /* Init the list of queues. */
        for ( k = 0 ; k < e->nr_queues ; k++ )
            queues[k] = &e->queues[k];
        naq = e->nr_queues - 1;
        queues[ myqid ] = queues[ naq ];
                        
        /* while i can still get a pair... */
        /* printf("runner_run: runner %i paSSEd barrier, getting pairs...\n",r->id); */
        while ( myq->next < myq->count || naq > 0 ) {
        
            /* Try to get a pair from my own queue. */
            TIMER_TIC
            if ( myq->next == myq->count || ( t = queue_get( myq , r->id , 0 ) ) == NULL ) {
            
                /* Clean up the list of queues. */
                count = myq->count - myq->next;
                for ( k = 0 ; k < naq ; k++ ) {
                    count += queues[k]->count - queues[k]->next;
                    if ( queues[k]->next == queues[k]->count )
                        queues[k--] = queues[--naq];
                    }
                        
                /* If there are no queues left, go back to go, do not collect 200 FLOPs. */
                if ( naq == 0 )
                    continue;
                    
                /* Otherwise, try to grab something from a random queue. */
                qid = rand_r( &myseed ) % naq;
                if ( ( t = queue_get( queues[qid] , r->id , 1 ) ) != NULL ) {
                
                    /* Add this task to my own queue. */
                    if ( !queue_insert( myq , t ) )
                        queue_insert( queues[qid] , t );
                
                    }
            
                /* If there are more queues than tasks, fall on sword. */
                if ( t == NULL && count <= r->id )
                    break;
                    
                }
                
            /* If I didn't get a task, try again, locking... */
            if ( t == NULL ) {
                
                /* Lock the mutex. */
                if ( pthread_mutex_lock( &s->tasks_mutex ) != 0 )
                    return error(runner_err_pthread);
                    
                /* Try again to get a pair... */
                if ( myq->next == myq->count || ( t = queue_get( myq , r->id , 0 ) ) == NULL ) {
                    count = myq->count - myq->next;
                    for ( k = 0 ; k < naq ; k++ ) {
                        count += queues[k]->count - queues[k]->next;
                        if ( queues[k]->next == queues[k]->count )
                            queues[k--] = queues[--naq];
                        }
                    if ( naq != 0 ) {
                        qid = rand_r( &myseed ) % naq;
                        if ( ( t = queue_get( queues[qid] , r->id , 1 ) ) != NULL ) {
                            if ( !queue_insert( myq , t ) )
                                queue_insert( queues[qid] , t );
                            }
                        }
                    }
                    
                /* If no pair, wait... */
                if ( count > 0 && t == NULL )    
                    if ( pthread_cond_wait( &s->tasks_avail , &s->tasks_mutex ) != 0 )
                        return error(runner_err_pthread);
                        
                /* Unlock the mutex. */
                if ( pthread_mutex_unlock( &s->tasks_mutex ) != 0 )
                    return error(runner_err_pthread);
                
                /* Skip back to the top of the queue if empty-handed. */
                if ( t == NULL )
                    continue;
                
                }
            TIMER_TOC(runner_timer_queue);
                
            /* Check task type... */
            switch ( t->type ) {
                case task_type_sort:
                    TIMER_TIC_ND
                    if ( s->verlet_rebuild && !( e->flags & engine_flag_unsorted ) )
                        if ( runner_dosort( r , &s->cells[ t->i ] , t->flags ) < 0 )
                            return error(runner_err);
                    s->cells_taboo[ t->i ] = 0;
                    TIMER_TOC(runner_timer_sort);
                    break;
                case task_type_self:
                    TIMER_TIC_ND
                    if ( runner_doself( r , &s->cells[ t->i ] ) < 0 )
                        return error(runner_err);
                    s->cells_taboo[ t->i ] = 0;
                    TIMER_TOC(runner_timer_self);
                    break;
                case task_type_pair:
                    TIMER_TIC_ND
                    if ( e->flags & engine_flag_unsorted ) {
                        if ( runner_dopair_unsorted( r , &s->cells[ t->i ] , &s->cells[ t->j ] ) < 0 )
                            return error(runner_err);
                        }
                    else {
                        if ( runner_dopair( r , &s->cells[ t->i ] , &s->cells[ t->j ] , t->flags ) < 0 )
                            return error(runner_err);
                        }
                    s->cells_taboo[ t->i ] = 0;
                    s->cells_taboo[ t->j ] = 0;
                    TIMER_TOC(runner_timer_pair);
                    break;
                default:
                    return error(runner_err_tasktype);
                }
                
            /* Unlock any dependent tasks. */
            for ( k = 0 ; k < t->nr_unlock ; k++ )
                __sync_fetch_and_sub( &t->unlock[k]->wait , 1 );

            /* Bing! */
            if ( pthread_mutex_lock( &s->tasks_mutex ) != 0 )
                return error(runner_err_pthread);
            if ( pthread_cond_broadcast( &s->tasks_avail ) != 0 )
                return error(runner_err_pthread);
            if ( pthread_mutex_unlock( &s->tasks_mutex ) != 0 )
                return error(runner_err_pthread);

            }

        /* give the reaction count */
        // printf("runner_run: last count was %u.\n",runner_rcount);
        r->err = acc;
            
        /* did things go wrong? */
        /* printf("runner_run: runner %i done pairs.\n",r->id); fflush(stdout); */
        if ( err < 0 )
            return error(runner_err_space);
    
        /* Bing! */
        if ( pthread_mutex_lock( &s->tasks_mutex ) != 0 )
            return error(runner_err_pthread);
        if ( pthread_cond_broadcast( &s->tasks_avail ) != 0 )
            return error(runner_err_pthread);
        if ( pthread_mutex_unlock( &s->tasks_mutex ) != 0 )
            return error(runner_err_pthread);

        }

    }

    
/**
 * @brief Initialize the runner associated to the given engine and
 *      attach it to an SPU.
 * 
 * @param r The #runner to be initialized.
 * @param e The #engine with which it is associated.
 * @param id The ID of this #runner.
 * 
 * @return #runner_err_ok or < 0 on error (see #runner_err).
 *
 * If @c CELL is not defined, this routine will fail!
 */

int runner_init_SPU ( struct runner *r , struct engine *e , int id ) {

#ifdef CELL
    static void *data = NULL;
    static int size_data = 0;
    void *finger;
    int nr_pots = 0, size_pots = 0, *pots, i, j, k, l;
    struct potential *p;
    unsigned int buff;

    /* make sure the inputs are ok */
    if ( r == NULL || e == NULL )
        return error(runner_err_null);
        
    /* remember who i'm working for */
    r->e = e;
    r->id = id;
    
    /* if this has not been done before, init the runner data */
    if ( data == NULL ) {

        /* run through the potentials and count them and their size */
        for ( i = 0 ; i < e->max_type ; i++ )
            for ( j = i ; j < e->max_type ; j++ )
                if ( e->p[ i * e->max_type + j ] != NULL ) {
                    nr_pots += 1;
                    size_pots += e->p[ i * e->max_type + j ]->n + 1;
                    }

        /* the main data consists of a pointer to the cell data (64 bit),
           the nr of cells (int), the cutoff (float), the width of
           each cell (float[3]), the max nr of types (int)
           and an array of size max_type*max_type of offsets (int) */
        size_data = sizeof(void *) + sizeof(int) + 4 * sizeof(float) + sizeof(int) * ( 1 + e->max_type*e->max_type );

        /* stretch this data until we are aligned to 8 bytes */
        size_data = ( size_data + 7 ) & ~7;

        /* we then append nr_pots potentials consisting of three floats (alphas) */
        /* and two ints (n and flags) */
        size_data += nr_pots * ( 3 * sizeof(float) + 2 * sizeof(int) );

        /* finally, we append the data of each interval of each potential */
        /* which consists of eight floats */
        size_data += size_pots * sizeof(float) * potential_chunk;

        /* raise to multiple of 128 */
        size_data = ( size_data + 127 ) & ~127;

        /* allocate memory for the SPU data */
        if ( ( data = malloc_align( size_data , 7 ) ) == NULL )
            return error(runner_err_malloc);

        /* fill-in the engine data (without the pots) */
        finger = data;
        *((unsigned long long *)finger) = 0; finger += sizeof(unsigned long long);
        *((int *)finger) = e->s.nr_cells; finger += sizeof(int);
        *((float *)finger) = e->s.cutoff; finger += sizeof(float);
        *((float *)finger) = e->s.h[0]; finger += sizeof(float);
        *((float *)finger) = e->s.h[1]; finger += sizeof(float);
        *((float *)finger) = e->s.h[2]; finger += sizeof(float);
        *((int *)finger) = e->max_type; finger += sizeof(int);
        pots = (int *)finger; finger += e->max_type * e->max_type * sizeof(int);
        for ( i = 0 ; i < e->max_type*e->max_type ; i++ )
            pots[i] = 0;

        /* move the finger until we are at an 8-byte boundary */
        finger = (void *)( ( (unsigned long long)finger + 7 ) & ~7 );

        /* loop over the potentials */
        for ( i = 0 ; i < e->max_type ; i++ )
            for ( j = i ; j < e->max_type ; j++ )
                if ( pots[ i * e->max_type + j ] == 0 && e->p[ i * e->max_type + j ] != NULL ) {
                    p = e->p[ i * e->max_type + j ];
                    for ( k = 0 ; k < e->max_type*e->max_type ; k++ )
                        if ( e->p[k] == p )
                            pots[k] = finger - data;
                    *((int *)finger) = p->n; finger += sizeof(int);
                    *((int *)finger) = p->flags; finger += sizeof(int);
                    *((float *)finger) = p->alpha[0]; finger += sizeof(float);
                    *((float *)finger) = p->alpha[1]; finger += sizeof(float);
                    *((float *)finger) = p->alpha[2]; finger += sizeof(float);
                    /* loop explicitly in case FPTYPE is not float. */
                    for ( k = 0 ; k <= p->n ; k++ ) {
                        for ( l = 0 ; l < potential_chunk ; l++ ) {
                            *((float *)finger) = p->c[k*potential_chunk+l];
                            finger += sizeof(float);
                            }
                        }
                    }

        /* raise to multiple of 128 */
        finger = (void *)( ( (unsigned long long)finger + 127 ) & ~127 );

        /* if the effective size is smaller than the allocated size */
        /* (e.g. duplicate potentials), be clean and re-allocate the data */
        if ( finger - data < size_data ) {
            size_data = finger - data;
            if ( ( data = realloc_align( data , size_data , 7 ) ) == NULL )
                return error(runner_err_malloc);
            }

        /* say something about it all */
        /* printf("runner_init: initialized data with %i bytes.\n",size_data); */

        } /* init runner data */

    /* remember where the data is */
    r->data = data;

    /* allocate and set the cell data */
    if ( ( r->celldata = (struct celldata *)malloc_align( ceil128( sizeof(struct celldata) * r->e->s.nr_cells ) , 7 ) ) == NULL )
        return error(runner_err_malloc);
    *((unsigned long long *)data) = (unsigned long long)r->celldata;

    /* get a handle on an SPU */
    if ( ( r->spe = spe_context_create(0, NULL) ) == NULL )
        return error(runner_err_spe);

    /* load the image onto the SPU */
    if ( spe_program_load( r->spe , &runner_spu ) != 0 )
        return error(runner_err_spe);

    /* dummy function that just starts the SPU... */
    int dummy ( struct runner *r ) {
        return spe_context_run( r->spe , &(r->entry) , 0 , r->data , (void *)(unsigned long long)size_data , NULL );
        }

    /* start the SPU with a pointer to the data */
    r->entry = SPE_DEFAULT_ENTRY;
	if (pthread_create(&r->spe_thread,NULL,(void *(*)(void *))dummy,r) != 0)
		return error(runner_err_pthread);

    /* wait until the SPU is ready... */
    if ( spe_out_intr_mbox_read( r->spe , &buff , 1 , SPE_MBOX_ALL_BLOCKING ) < 1 )
        return runner_err_spe;

    /* start the runner. */
    if ( e->flags & engine_flag_tuples ) {
	    if (pthread_create(&r->thread,NULL,(void *(*)(void *))runner_run_cell_tuples,r) != 0)
		    return error(runner_err_pthread);
        }
    else {
	    if (pthread_create(&r->thread,NULL,(void *(*)(void *))runner_run_cell,r) != 0)
		    return error(runner_err_pthread);
        }

    /* all is well... */
    return runner_err_ok;
    
#else
        
    /* if not compiled for cell, then this option is not available. */
    return error(runner_err_unavail);
        
#endif
        
    }
    
    
/**
 * @brief Initialize the runner associated to the given engine.
 * 
 * @param r The #runner to be initialized.
 * @param e The #engine with which it is associated.
 * @param id The ID of this #runner.
 * 
 * @return #runner_err_ok or < 0 on error (see #runner_err).
 */

int runner_init ( struct runner *r , struct engine *e , int id ) {

    #if defined(HAVE_SETAFFINITY) && !defined(CELL)
        cpu_set_t cpuset;
    #endif

    /* make sure the inputs are ok */
    if ( r == NULL || e == NULL )
        return error(runner_err_null);
        
    /* remember who i'm working for */
    r->e = e;
    r->id = id;
    
    /* init the thread using tasks. */
	if ( pthread_create( &r->thread , NULL , (void *(*)(void *))runner_run , r ) != 0 )
		return error(runner_err_pthread);
    
    /* If we can, try to restrict this runner to a single CPU. */
    #if defined(HAVE_SETAFFINITY) && !defined(CELL)
        if ( e->flags & engine_flag_affinity ) {
        
            /* Set the cpu mask to zero | r->id. */
            CPU_ZERO( &cpuset );
            CPU_SET( r->id , &cpuset );

            /* Apply this mask to the runner's pthread. */
            if ( pthread_setaffinity_np( r->thread , sizeof(cpu_set_t) , &cpuset ) != 0 )
                return error(runner_err_pthread);

            }
    #endif
    
    /* all is well... */
    return runner_err_ok;
    
    }
