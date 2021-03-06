
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <unistd.h>
#include <errno.h>
#include <memory.h>
#include <stdarg.h>
#include <time.h>
#include <ctype.h>
#include <math.h>
#include <pwd.h>

#include <libpq-fe.h>
#include "database.h"


//====================================================
//	UTILS
//====================================================

int getldate(time_t tt)
{
	struct tm *pnow = localtime( &tt ); /* Convert to local time. */
	
	return pnow->tm_mday + 100*(pnow->tm_mon+1) + 10000*(pnow->tm_year+1900);
}

int getlnow()
{	
	return getldate(time(NULL));
}

unsigned long long getCurrentUTS()
{
	unsigned long long uts = 0ULL;
	unsigned long long uty = 0ULL;
	time_t tt;
	
	time( &tt );                /* Get time as long integer. */
	struct tm *pnow = localtime( &tt ); /* Convert to local time. */
	
	uts = pnow->tm_sec +100*pnow->tm_min + 10000*pnow->tm_hour + 1000000*pnow->tm_mday + 100000000L*(pnow->tm_mon+1);
	
	uty = 10000000000ULL*(pnow->tm_year+1900);
	uts += uty;
	
	return uts;
}


char *formatTS(long long ts, char *buffer)
{
	int time = ts % 1000000;
	int seconds = time % 100;
	int restant = time / 100;
	int minutes = restant % 100;
	int hour = restant / 100;
	int date = ts / 1000000;
	int day = date % 100;
	restant = date / 100;
	int month = restant % 100;
	int year = restant / 100;

	sprintf(buffer, "%4d-%02d-%02d  %02d:%02d:%02d", year, month, day, hour, minutes, seconds);

	return buffer;
}

//====================================================
//	POSTGRES
//====================================================

void pgCloseConn(PGconn *conn)
{
	PQfinish(conn);
}

PGconn *pgOpenConn(const char *pdbname, const char *pdbusername, const char *pdbpassword, char *buffer_error)
{
PGconn *pgConn = NULL;
char connect_string[256];

	sprintf(connect_string, "dbname=%s user=%s password=%s",
		pdbname,
		pdbusername,
		pdbpassword);
		
	pgConn = PQconnectdb(connect_string);

	/*
    * check to see that the backend connection was successfully made
    */
	if (PQstatus(pgConn) == CONNECTION_BAD)
	{
		if (buffer_error != NULL)
		{
			strcpy(buffer_error, PQerrorMessage(pgConn));
		}
		return NULL;
	}
	else
	{
		return pgConn;
	}
}

bool pgExec(PGconn *pgConn, int *pnbc, const char *pCommand)
{
bool bOK = false;
PGresult   *pgres = NULL;

	*pnbc = 0;
        pgres = PQexec(pgConn, pCommand);
        if (pgres == NULL || PQresultStatus(pgres) != PGRES_COMMAND_OK)
        {
		printf("\nSQL_ERROR: '%s'  [%s]\n", pCommand, PQerrorMessage(pgConn));
		if (pgres != NULL) PQclear(pgres);
		bOK = false;
        }
	else
	{
		*pnbc = atoi(PQcmdTuples(pgres));
		bOK = true;
	}

	return bOK;
}

PGresult *pgQuery(PGconn *pgConn, const char *pQuery)
{
	PGresult   *pgres = NULL;

	pgres = PQexec(pgConn, pQuery);
	if (pgres == NULL)
	{
		printf("pgres is NULL\n");
	}
	else
	{
		ExecStatusType est = PQresultStatus(pgres);
		if (est != PGRES_TUPLES_OK)
		{
			printf("SQL_ERROR: '%s'  [%s]\n", pQuery, PQerrorMessage(pgConn));
			if (pgres != NULL) PQclear(pgres);
			pgres = NULL;
		}
	}
	return pgres;
}

bool pgExecFormat(PGconn *pgConn, int *pnbc, char *szFormat, ...)
{
	char buffer[4800];
	va_list args;

	va_start( args, szFormat );
	vsprintf( buffer, szFormat, args);
	va_end( args );
	*pnbc = 0;
	
	return pgExec(pgConn, pnbc, buffer);
}

bool pgBeginTransaction(PGconn *pgConn)
{
	int nbc = 0;
	
	return pgExec(pgConn, &nbc, "BEGIN");
}

bool pgCommitTransaction(PGconn *pgConn)
{
	int nbc = 0;
	
	return pgExec(pgConn, &nbc, "COMMIT");
}

bool pgRollbackTransaction(PGconn *pgConn)
{
	int nbc = 0;
	
	return pgExec(pgConn, &nbc, "ROLLBACK");
}

long long pgGetSequo(PGconn *pgConn, const char *psequo)
{
char query[256];
long long lval = -1;
PGresult   *pgres = NULL;

	sprintf(query, "select nextval('%s')", psequo);

	pgres = PQexec(pgConn, query);
	if (pgres != NULL)
	{
		if (PQresultStatus(pgres) != PGRES_TUPLES_OK)
		{
			printf("Database ERROR: %s\n", PQerrorMessage(pgConn));
			lval = -2;
		}
		else
		{
			if (PQntuples(pgres) > 0)
			{
				lval = atoll(PQgetvalue(pgres, 0, 0));
			}
		}
		PQclear(pgres);
	}
	else
	{
		printf("pgGetSequo(): NULL PGresult returned.\n");
		lval = -3;
	}

	return lval;
}

int pgGetCount(PGconn *pgConn, const char *ptable, const char *pwhere)
{
	char query[512];
	int count = -1;
	PGresult   *pgres = NULL;

	if (pwhere != NULL && strlen(pwhere) > 0)
		sprintf(query, "select count(*) from %s where %s", ptable, pwhere);
	else
		sprintf(query, "select count(*) from %s", ptable);

	pgres = PQexec(pgConn, query);
	if (pgres != NULL)
	{
		if (PQresultStatus(pgres) != PGRES_TUPLES_OK)
		{
			printf("Database ERROR: %s\n", PQerrorMessage(pgConn));
		}
		else
		{
			if (PQntuples(pgres) > 0)
			{
				count = atoi(PQgetvalue(pgres, 0, 0));
			}
		}
		PQclear(pgres);
	}
	else
	{
		printf("pgGetCount(): NULL PGresult returned.\n");
	}

	return count;
}


double pgGetDouble(PGconn *pgConn, const char *pQuery)
{
	double dValue = 0.0;
	PGresult   *pgres = NULL;


	pgres = PQexec(pgConn, pQuery);
	if (pgres != NULL)
	{
		if (PQresultStatus(pgres) != PGRES_TUPLES_OK)
		{
			printf("Database ERROR: %s\n", PQerrorMessage(pgConn));
		}
		else
		{
			if (PQntuples(pgres) > 0)
			{
				dValue = atoi(PQgetvalue(pgres, 0, 0));
			}
		}
		PQclear(pgres);
	}
	else
	{
		printf("pgGetDouble(): NULL PGresult returned.\n");
	}

	return dValue;
}

int pgGetMax(PGconn *pgConn, const char *pColumn, const char *ptable, const char *pwhere)
{
	char query[256];
	int max = -1;
	PGresult   *pgres = NULL;

	if (pwhere != NULL && strlen(pwhere) > 0)
		sprintf(query, "select max(%s) from %s where %s", pColumn, ptable, pwhere);
	else
		sprintf(query, "select max(%s) from %s", pColumn, ptable);

	pgres = PQexec(pgConn, query);
	if (pgres != NULL)
	{
		if (PQresultStatus(pgres) != PGRES_TUPLES_OK)
		{
			printf("Database ERROR: %s\n", PQerrorMessage(pgConn));
		}
		else
		{
			if (PQntuples(pgres) > 0)
			{
				max = atoi(PQgetvalue(pgres, 0, 0));
			}
		}
		PQclear(pgres);
	}
	else
	{
		printf("pgGetMax(): NULL PGresult returned.\n");
	}

	return max;
}

double pgGetSum(PGconn *pgConn, const char *pColumn, const char *ptable, const char *pwhere)
{
	char query[256];
	double sum = 0.0;
	PGresult   *pgres = NULL;

	if (pwhere != NULL && strlen(pwhere) > 0)
		sprintf(query, "select sum(%s) from %s where %s", pColumn, ptable, pwhere);
	else
		sprintf(query, "select sum(%s) from %s", pColumn, ptable);

	pgres = PQexec(pgConn, query);
	if (pgres != NULL)
	{
		if (PQresultStatus(pgres) != PGRES_TUPLES_OK)
		{
			printf("Database ERROR: %s\n", PQerrorMessage(pgConn));
		}
		else
		{
			if (PQntuples(pgres) > 0)
			{
				sum = atof(PQgetvalue(pgres, 0, 0));
			}
		}
		PQclear(pgres);
	}
	else
	{
		printf("pgGetSum(): NULL PGresult returned.\n");
	}

	return sum;
}

//====================================================
//	SIMULATION
//====================================================

int insertSimulation(PGconn *pgConn, int epochs, int frames,
                     double gravity, double rebound, double transfer, double factor, double cf1, double cf2,
                     int mmin, int mmax, int vmin, int vmax, int bc, int pc, int score)
{
	int nbc = 0, sid = -1;
	unsigned long long uts = getCurrentUTS();

	if (pgBeginTransaction(pgConn))
	{
		sid = pgGetSequo(pgConn, "g4s.simulation_sequo");
		if (sid < 0)
		{
			pgRollbackTransaction(pgConn);
		}
		else
		{
			pgExecFormat(pgConn, &nbc, "insert into g4s.simulation values (%d, %llu, %d, %d, %6.2f, %6.2f, %6.2f, %6.2f, %6.2f, %6.2f, %d, %d, %d, %d, %d, %d, %d, 0, 0)",
sid, uts, epochs, frames, gravity, rebound, transfer, factor, cf1, cf2, mmin, mmax, vmin, vmax, bc, pc, score);
			if (nbc == 1)
				pgCommitTransaction(pgConn);
			else
			{
				pgRollbackTransaction(pgConn);
				sid = -1;
			}
		}
	}
	return sid;
}

bool insertFrame(PGconn *pgConn, int sid, int epoch, int pid, int x, int y, double radius, const char *color)
{
int nbc = 0;
	pgExecFormat(pgConn, &nbc, "insert into g4s.frame values (%d, %d, %d, %d, %d, %6.2f, '%s')",
		sid, epoch, pid, x, y, radius, color);
	return nbc == 1;
}

