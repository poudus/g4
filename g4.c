
#include <string.h>
#include <stdio.h>
#include <ctype.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <time.h>

typedef struct
{
int     epochs, frames, planets, ec0, ecf, bc, pc;
double  gravity, rebound;
} SIMULATION;

typedef struct
{
int     id, mass, red, green, blue;
} PLANET;

typedef struct
{
int     pid, epoch, np;
double  x, y, vx, vy;
} FRAME;

//
//	GLOBALS
//

#define NB_MAX_FRAMES	65536

PLANET	planet[16];
FRAME	frame[NB_MAX_FRAMES];
int     red[6], green[6], blue[6];

//--------------------------------

PLANET* init_planet(int id, int mass, int red, int green, int blue, PLANET *planet)
{
	planet->id = id;
	planet->mass = mass;
	planet->red = red;
	planet->green = green;
	planet->blue = blue;
	return planet;
}

void init_planets(int planets, int m1, int m2, PLANET *planet)
{

	for (int p = 0 ; p < planets ; p++)
	{
		init_planet(p, m1 + rand() % (m2 - m1 + 1), red[p], green[p], blue[p], planet);
		planet++;
	}
}

FRAME* init_frame(int pid, int epoch, double x, double y, double vx, double vy, FRAME* frame)
{
	frame->pid = pid;
	frame->np = 0;
	frame->epoch = epoch;
	frame->x = x;
	frame->y = y;
	frame->vx = vx;
	frame->vy = vy;
	return frame;
}

void init_frames(int frames, int offset, int size, int v1, int v2, FRAME* frame)
{
	for (int f = 0 ; f < frames ; f++)
	{
		init_frame(f, 0, offset + rand() % (size-2*offset), offset + rand() % (size-2*offset),
			v1 + rand() % (v2 - v1 + 1), v1 + rand() % (v2 - v1 + 1), frame);
		frame++;
	}
}

double ec(int frames, FRAME* frame)
{
double ec = 0.0;

	for (int f = 0 ; f < frames ; f++)
		ec += planet[frame[f].pid].mass * (frame[f].vx * frame[f].vx + frame[f].vy * frame[f].vy);
	return ec;
}

void	attraction(int pid, int x, int y, int frames, FRAME* frame, double *mgx, double* mgy)
{
	*mgx = *mgy = 0.0;

	for (int f = 0 ; f < frames ; f++)
	{
		if (frame[f].pid != pid)
		{
			*mgx += planet[frame[f].pid].mass / (x-frame[f].x)*(x-frame[f].x);
			*mgy += planet[frame[f].pid].mass / (y-frame[f].y)*(y-frame[f].y);
		}
	}
}

double distance(FRAME* f1, FRAME* f2)
{
	return sqrt((f2->x - f1->x) * (f2->x - f1->x) + (f2->y - f1->y) * (f2->y - f1->y));
}

double angle(int x1, int y1, int x2, int y2, int vx1, int vy1, double* vr)
{
double ar = atan2(y2 - y1, x2 - x1);
double av = atan2(vy1, vx1);

	*vr = cos(ar-av) * sqrt(vx1*vx1 + vy1*vy1);

	return ar-av;
}

void bc(int m1, int m2, double d, FRAME* f1, FRAME* f2, double rebound)
{
    double vr1, vr2;
    
    double a1 = angle(f1->x, f1->y, f2->x, f2->y, f1->vx, f1->vy, &vr1);
    double a2 = angle(f2->x, f2->y, f1->x, f1->y, f2->vx, f2->vy, &vr2);

//printf("PC %d/%4d - %d/%4d   a1 = %5.2f  a2 = %5.2f  [ %5.2f  %5.2f ]\n",
//        f1->pid, m1, f2->pid, m2, 180.0*a1/3.14159, 180.0*a2/3.14159, vr1, vr2);

//printf("PC x1 = %5.2f  y1 = %5.2f   x2 = %5.2f  y2 = %5.2f\n", f1->x, f1->y, f2->x, f2->y);

	//g1 = rebound * qvr2 / (100.0 * (m1+m2));
	double g1x = (vr1 * (f2->x - f1->x)) / d;
	double g1y = (vr1 * (f2->y - f1->y)) / d;
	//g2 = rebound * m1 * vr1 / (100.0 * (m1+m2));
	double g2x = (vr2 * (f1->x - f2->x)) / d;
	double g2y = (vr2 * (f1->y - f2->y)) / d;

//printf("PC g1x = %5.2f  g1y = %5.2f   g2x = %5.2f  g2y = %5.2f\n", g1x, g1y, g2x, g2y);

//printf("before v1x = %6.2f  v1y = %6.2f   v2x = %6.2f  v2y = %6.2f\n", f1->vx, f1->vy, f2->vx, f2->vy);

    f1->vx = (m1 * f1->vx + 0.01 * rebound * m2 * g2x) / (m1 + m2);
	f1->vy = (m1 * f1->vy + 0.01 * rebound * m2 * g2y) / (m1 + m2);

	f2->vx = (m2 * f2->vx + 0.01 * rebound * m1 * g1x) / (m1 + m2);
	f2->vy = (m2 * f2->vy + 0.01 * rebound * m1 * g1y) / (m1 + m2);
    
//printf("after  v1x = %6.2f  v1y = %6.2f   v2x = %6.2f  v2y = %6.2f\n", f1->vx, f1->vy, f2->vx, f2->vy);

    f1->x += f1->vx;
	f1->y += f1->vy;

	f2->x += f2->vx;
	f2->y += f2->vy;
}

typedef struct
{
    int m1, m2, np;
    double x1, y1, vx1, vy1, x2, y2, vx2, vy2;
} COALESCE;

bool run_simulation(int planets, int size, int mmin, int mmax, int vmin, int vmax, double gravity, double rebound,
	SIMULATION* simulation)
{
int	offset = 100, first_frame, last_frame, nf, m1, m2, p1, p2;
double	mgx, mgy, radius, r1, r2, d;
bool	stop = false;
COALESCE    coalesce[8];

	simulation->planets = planets;
	simulation->gravity = gravity;
	simulation->rebound = rebound;
	simulation->bc = 0;
	simulation->pc = 0;

	init_planets(planets, mmin, mmax, &planet[0]);
	init_frames(planets, offset, size, vmin, vmax, &frame[0]);

	simulation->epochs = 1;
	simulation->frames = planets;

	simulation->ec0 = ec(planets, &frame[0]);
	printf("%d planets, gravity = %5.2f, rebound = %5.2f, mass = [ %d - %d ], v0 = [ %d - %d ], ec0 = %d K\n",
		planets, gravity, rebound, mmin, mmax, vmin, vmax, simulation->ec0/1000);

	first_frame = 0;
	last_frame = planets;
    int idx_new_planet = planets;

	while (!stop)
	{
		// gravity

		for (int f = first_frame ; f < last_frame ; f++)
		{
			attraction(frame[f].pid, frame[f].x, frame[f].y,
				last_frame-first_frame, &frame[first_frame], &mgx, &mgy);

			nf = simulation->frames;

			frame[nf].vx = frame[f].vx + (gravity*mgx)/(planet[frame[f].pid].mass*100.0);
			frame[nf].vy = frame[f].vy + (gravity*mgy)/(planet[frame[f].pid].mass*100.0);

			frame[nf].x = frame[f].x + frame[nf].vx;
			frame[nf].y = frame[f].y + frame[nf].vy;


			frame[nf].epoch = simulation->epochs;
			frame[nf].pid = frame[f].pid;
			frame[nf].np = 0;

			simulation->frames++;
		}

		first_frame = last_frame;
		last_frame = simulation->frames;
		simulation->epochs++;

		// border collisions

		for (int f = first_frame ; f < last_frame ; f++)
		{
			radius = sqrt(planet[frame[f].pid].mass);
			if (frame[f].x < radius) // LEFT
			{
				frame[f].x = 2*radius - frame[f].x;
				frame[f].vx = -frame[f].vx;
				simulation->bc++;
			}
			else if (frame[f].x > (size-radius)) // RIGHT
			{
				frame[f].x = 2*(size-radius) - frame[f].x;
				frame[f].vx = -frame[f].vx;
				simulation->bc++;
			}
			else if (frame[f].y < radius) // TOP
			{
				frame[f].y = 2*radius - frame[f].y;
				frame[f].vy = -frame[f].vy;
				simulation->bc++;
			}
			else if (frame[f].y > (size-radius)) // BOTTOM
			{
				frame[f].y = 2*(size-radius) - frame[f].y;
				frame[f].vy = -frame[f].vy;
				simulation->bc++;
			}
		}

		// planet collisions
        int nb_coalesce = 0;
		for (int f1 = first_frame ; f1 < last_frame ; f1++)
		{
			p1 = frame[f1].pid;
			m1 = planet[frame[f1].pid].mass;
			r1 = sqrt(m1);
			for (int f2 = first_frame ; f2 < last_frame ; f2++)
			{
				p2 = frame[f2].pid;
				m2 = planet[frame[f2].pid].mass;
				r2 = sqrt(m2);
				if (f1 < f2)
				{
					d = distance(&frame[f1], &frame[f2]);
					if (d < r1 + r2)
					{
printf("COLLISION    e = %4d  %4d/%d and %4d/%d\n",
       simulation->epochs, f1, frame[f1].pid, f2, frame[f2].pid);

                        bc(m1, m2, d, &frame[f1], &frame[f2], rebound);
                        d = distance(&frame[f1], &frame[f2]);
                        if (d < r1 + r2) // coalescence
                        {
                            coalesce[nb_coalesce].m1 = m1;
                            coalesce[nb_coalesce].m2 = m2;
                            coalesce[nb_coalesce].np = idx_new_planet;
                            coalesce[nb_coalesce].x1 = frame[f1].x;
                            coalesce[nb_coalesce].y1 = frame[f1].y;
                            coalesce[nb_coalesce].vx1 = frame[f1].vx;
                            coalesce[nb_coalesce].vy1 = frame[f1].vy;
                            coalesce[nb_coalesce].x2 = frame[f2].x;
                            coalesce[nb_coalesce].y2 = frame[f2].y;
                            coalesce[nb_coalesce].vx2 = frame[f2].vx;
                            coalesce[nb_coalesce].vy2 = frame[f2].vy;
                            nb_coalesce++;
                            planet[idx_new_planet].id = idx_new_planet;
                            planet[idx_new_planet].mass = m1 + m2;
                            planet[idx_new_planet].red = (m1*planet[p1].red+m2*planet[p2].red)/(m1+m2);
                            planet[idx_new_planet].green = (m1*planet[p1].green+m2*planet[p2].green)/(m1+m2);
                            planet[idx_new_planet].blue = (m1*planet[p1].blue+m2*planet[p2].blue)/(m1+m2);
                            frame[f1].np = frame[f2].np = idx_new_planet;
printf("COALESCENCE  e = %4d  %4d/%d and %4d/%d   np = %d\n",
       simulation->epochs, f1, frame[f1].pid, f2, frame[f2].pid, idx_new_planet);
                            idx_new_planet++;
                        }
						simulation->pc++;
					}
				}
			}
		}
		if (nb_coalesce > 0) // coalesce
        {
            FRAME   tmpf[16];
            int     nbtmpf = 0;
            
            for (int f = first_frame ; f < last_frame ; f++)
            {
                if (frame[f].np == 0)
                {
                    memcpy(&tmpf[nbtmpf], &frame[f], sizeof(FRAME));
                    nbtmpf++;
                }
            }
            for (int tf = 0 ; tf < nbtmpf ; tf++)
                memcpy(&frame[first_frame+tf], &tmpf[tf], sizeof(FRAME));
            last_frame = first_frame + nbtmpf;
            for (int cf = 0 ; cf < nb_coalesce ; cf++)
            {
                frame[last_frame+cf].epoch = simulation->epochs;
                frame[last_frame+cf].pid = coalesce[cf].np;
                frame[last_frame+cf].np = 0;
                
                frame[last_frame+cf].x = (coalesce[cf].m1 * coalesce[cf].x1 + coalesce[cf].m2 * coalesce[cf].x2) / (coalesce[cf].m1 + coalesce[cf].m2);
                frame[last_frame+cf].y = (coalesce[cf].m1 * coalesce[cf].y1 + coalesce[cf].m2 * coalesce[cf].y2) / (coalesce[cf].m1 + coalesce[cf].m2);
                
                frame[last_frame+cf].vx = (coalesce[cf].m1 * coalesce[cf].vx1 + coalesce[cf].m2 * coalesce[cf].vx2) / (coalesce[cf].m1 + coalesce[cf].m2);
                frame[last_frame+cf].vy = (coalesce[cf].m1 * coalesce[cf].vy1 + coalesce[cf].m2 * coalesce[cf].vy2) / (coalesce[cf].m1 + coalesce[cf].m2);
                
                last_frame++;
            }
            planets = nbtmpf + nb_coalesce;
            printf("+++ planets = %d   nbtmpf = %d  nb_coalesce = %d  first frame = %d  last_frame = %d  nb frames = %d +++\n", planets, nbtmpf, nb_coalesce, first_frame, last_frame, simulation->frames);
            for (int f = first_frame ; f < last_frame ; f++)
            {
                printf("f%d : %d\n", f, frame[f].pid);
            }
            simulation->frames = last_frame;
        }
		// end conditions
		if (planets == 1 ||
            simulation->frames >= NB_MAX_FRAMES || simulation->epochs > 10000 || idx_new_planet > 15)
		{
			stop = true;
		}
	}

	return true;
}

//
// main
//
int main(int argc, char* argv[])
{
int	planets = 3, size = 800, m1 = 100, m2 = 900, v1 = 10, v2 = 30;
double	ec0 = 0.0, gravity = 5.0, rebound = 80.0;
SIMULATION simulation;

	red[0] = 255; red[1] = 0; red[2] = 0; red[3] = 255; red[4] = 255; red[5] = 0;
	green[0] = 0; green[1] = 255; green[2] = 0; green[3] = 255; green[4] = 0; green[5] = 255;
	blue[0] = 0; blue[1] = 0; blue[2] = 255; green[3] = 0; green[4] = 255; green[5] = 255;

	srand(time(NULL));

	if (argc > 1)
		planets = atoi(argv[1]);
	if (argc > 2)
		gravity = atof(argv[2]);
	if (argc > 3)
		rebound = atof(argv[3]);
	if (argc > 4)
		m2 = atoi(argv[4]);
	if (argc > 5)
		v2 = atoi(argv[5]);
/*
double vr1, vr2;
double a1 = angle(1, 3, 3, 2, 1, 1, &vr1);
double a2 = angle(3, 2, 1, 3, 0, 1, &vr2);
printf("  a1 = %5.2f    a2 = %5.2f\n", a1, a2);
printf("cos1 = %5.2f  cos2 = %5.2f\n", cos(a1), cos(a2));
printf(" vr1 = %5.2f   vr2 = %5.2f\n", vr1, vr2);
exit(0);
*/

FRAME f1, f2;

f1.x = 1;
f1.y = 1;
f1.vx = -1;
f1.vy = 1;
f2.x = 3;
f2.y = 1;
f2.vx = -1;
f2.vy = 0;

    //bc(1, 1, 2.0, &f1, &f2, 80.0);
    //exit(0);
                        
	if (run_simulation(planets, size, m1, m2, v1, v2, gravity, rebound, &simulation))
		printf("%d epochs, %d frames  (%d bc) (%d pc)\n",
			simulation.epochs, simulation.frames, simulation.bc, simulation.pc);
	else
		printf("ERROR\n");
}

