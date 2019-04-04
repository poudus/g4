
#include <string.h>
#include <stdio.h>
#include <ctype.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <time.h>

typedef struct
{
} SIMULATION;

typedef struct
{
int	id, mass, red, green, blue;
} PLANET;

typedef struct
{
int	sid, pid, epoch;
double	x, y, vx, vy;
} FRAME;

PLANET* init_planet(int id, int mass, int red, int green, int blue, PLANET *planet)
{
	planet->id = id;
	planet->mass = mass;
	planet->red = red;
	planet->green = green;
	planet->blue = blue;
	return planet;
}

FRAME* init_frame(int sid, int pid, int epoch, double x, double y, double vx, double vy, FRAME* frame)
{
	frame->sid = sid;
	frame->pid = pid;
	frame->epoch = epoch;
	frame->x = x;
	frame->y = y;
	frame->vx = vx;
	frame->vy = vy;
	return frame;
}

//
// main
//
int main(int argc, char* argv[])
{
int	planets = 3, v0max = 30, depth = 2;
double	ec0 = 0.0, gravity = 5.0, rebound = 80.0;

	if (argc > 1)
		planets = atoi(argv[1]);
	if (argc > 2)
		v0max = atoi(argv[2]);

	printf("%d planets, gravity = %5.2f, rebound = %5.2f, ec0 = %d\n", planets, gravity, rebound, ec0);

	{
		sprintf(database_name, "k%s", argv[1]);
		strcpy(buffer_error, argv[1]);
		buffer_error[2] = 0;
		size = atoi(buffer_error);
		//printf("\nsize         = %d\n", size);
		pgConn = pgOpenConn(database_name, "k2", "", buffer_error);
		//printf("database connection     = %p %s\n", pgConn, database_name);
		//root[0] = 0;
		if (argc > 2)
			depth = atoi(argv[2]);
		if (argc > 3)
			min_count = atoi(argv[3]);
		if (argc > 4)
			min_elo_sum = atof(argv[4]);

		int nbi = 0, nbd = 0, len_root = strlen(root);
		bool bInsert = true; //(len_root == 0);
		//if (root[0] == '.') len_root = 0;
		int nb_book_moves = ComputeBookMoves(pgConn, &bmoves[0], size, min_count, depth, min_elo_sum);

		if (bInsert)
			printf("\n%d deleted\n", DeleteBookMoves(pgConn, depth));
		else
			printf("++\n");
		for (int bm = 0 ; bm < nb_book_moves ; bm++)
		{
			if (bInsert && InsertBookMove(pgConn, depth, &bmoves[bm])) nbi++;
			/*else if (len_root == 0 || strncmp(root, bmoves[bm].key, len_root) == 0)
			{
				printf("%s.%s   %6.2f %%   %5d  %4d - %-4d\n",
					bmoves[bm].key, bmoves[bm].move, bmoves[bm].ratio,
					bmoves[bm].count, bmoves[bm].win, bmoves[bm].loss);
				nbd++;
			}*/
		}
		if (bInsert)
			printf("%d/%d inserted\n\n", nbi, nb_book_moves);
		else
			printf("--\n%d book moves  (%d entries)\n\n", nbd, nb_book_moves);
	}
	else
	{
		printf("error.invalid-database-name\n");
		return -1;
	}
}

