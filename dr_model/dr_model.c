/* dr_model.c
 * Copyright (C) 2010 Battelle Memorial Institute
 * All Rights Reserved
 *
 * Author: David Chassin, Pacific Northwest National Laboratory
 * 
 * History:
 * 1.0 Created 21 May 2010 by David Chassin
 *
 * Notes: 
 * 1. This code has not yet been optimized and should not be 
 *    until an automated validation test is implemented
 * 2. The delta response is not yet implemented.
 *
 */

#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <memory.h>
#include "dr_model.h"

static char *errmsg = NULL;
static unsigned long nan[2]={0xffffffff, 0x7fffffff};

/* Get the last demand response model error 

	Returns:
		A pointer to the error message, if any, or NULL
*/
char *drm_error(void)
{
	return errmsg;
}

/* Create a demand response model

	Arguments:
		unsigned char L		the control band (see remarks)
		double eta			the initial demand (0.02 is typical)
		double phi			the initial duty cycle (0.5 is 50%)

	Returns:
		pointer to a DRMODEL structure or NULL on error

	Remarks:
		The time step of the model is 1/L.  This implies that it
		should take L timesteps for a device to go from the 'off'
		control bound at 0 to the 'on' control bound at L.

*/
DRMODEL *drm_create(unsigned short L, 
					double eta, 
					double phi)
{
	DRMODEL *drm = NULL;
	int x = 0;
	double r = 0;

	/* check the input values */
	if (eta<=-1 || eta>=1)
	{
		errmsg = "drm_create(): eta is invalid";
		return NULL;
	}
	if (phi<=0 || phi>=1)
	{
		errmsg = "drm_create(): phi is invalid";
		return NULL;
	}

	/* allocate memory for the model */
	drm = (DRMODEL*)malloc(sizeof(DRMODEL));
	if (drm==NULL)
	{
		errmsg = "drm_create(): memory allocation error";
		return NULL;
	}

	/* setup the buffers */
	drm->nbins = L;
	drm->off = (double*)malloc(sizeof(double)*L);
	drm->on = (double*)malloc(sizeof(double)*L);
	if (drm->off==NULL || drm->on==NULL)
	{
		errmsg = "drm_create(): memory allocation error";
		return NULL;
	}

	/* setup the initial population */
	r = phi/(1-phi);
	for (x=0; x<L; x++)
	{
		/* exponential distribution */
		if (eta==0)
		{
			drm->on[x] = eta * (1-phi) * exp(eta*(L-x+0.5)/r) / (exp(eta*L/r)-1);
			drm->off[x] = drm->on[x]/r;
		}

		/* uniform distribution */
		else
		{
			drm->on[x] = phi/L;
			drm->off[x] = (1-phi)/L;
		}
	}

	return drm;
}

/* Copy a demand response model 

	Arguments
		DRMODEL *from		the original to copy

	Returns
		A pointer to the new model 

	Remarks
		You must free the copy when you are done with it
*/
DRMODEL *drm_copy(DRMODEL *from)
{
	DRMODEL *copy = (DRMODEL*)malloc(sizeof(DRMODEL));
	if (copy==NULL)
	{
		errmsg = "drm_copy(): memory allocation error";
		return NULL;
	}
	copy->nbins = from->nbins;
	copy->on = (double*)malloc(from->nbins*sizeof(double));
	if (copy->on==NULL)
	{
		errmsg = "drm_copy(): memory allocation error";
		free(copy);
		return NULL;
	}
	copy->off = (double*)malloc(from->nbins*sizeof(double));
	if (copy->off==NULL)
	{
		errmsg = "drm_copy(): memory allocation error";
		free(copy->on);
		free(copy->off);
		return NULL;
	}
	memcpy(copy->on,from->on,from->nbins*sizeof(double));
	memcpy(copy->off,from->off,from->nbins*sizeof(double));
	return copy;
}

/* Make two models equal
	
	Arguments:
		DRMODEL *to			the model that is to be modified
		DRMODEL *from		the model that is to be copied

	Returns:
		1 on success, 0 on failure

	Remarks:
		This function only copies the population statistics.

*/
int drm_equal(DRMODEL *to, DRMODEL *from)
{
	if (to->nbins==from->nbins)
	{
		memcpy(to->on,from->on,from->nbins*sizeof(double));
		memcpy(to->off,from->off,from->nbins*sizeof(double));
		return 1;
	}
	else
	{
		errmsg = "drm_equal(): nbins don't match";
		return 0;
	}
}

/* Update a demand response model
  
	Arguments:
		DRMODEL *drm		the demand response model
		char delta			the band control shift
		double eta			the new demand (-1.0 to +1.0)
		double phi			the new duty cycle (0.0 to 1.0)

	Returns:
		The fraction of load that is on. QNAN if arguments are bad.

	Remarks:
		The time step of the model is invariant from one call to
		the next.  See drm_create for details.

		The model is updated so that if you intend to explore
		alternative, you must keep a copy of the model structure
		and reuse the original.  The version you want to commit
		is the only one you would keep for future updates.

*/
double drm_update(DRMODEL *drm, 
				  short delta, 
				  double eta, 
				  double phi)
{
	const double QNAN = *(double*)nan;
	double dnon[256], dnoff[256];
	double non=0, noff=0;
	int x;
	const int L = drm->nbins;
	double r;

	/* check the input values */
	if (delta!=0)
	{
		errmsg = "drm_update(): non-zero delta not supported yet";
		return QNAN;
	}
	if (eta<-1 || eta>1)
	{
		errmsg = "drm_update(): eta is invalid";
		return QNAN;
	}
	if (phi<=0 || phi>=1)
	{
		errmsg = "drm_update(): phi is invalid";
		return QNAN;
	}

	/* compute state changes */
	r = phi/(1-phi);
	for (x=0; x<L; x++)
	{
		/* normal mode */
		if (eta>=0)
		{
			/* on regime */
			if (x==L-1) /* upper boundary */
				dnon[x] = -drm->on[x] + (1-eta)*r*drm->off[x] + eta*drm->off[x];
			else /* on queue */
				dnon[x] = -drm->on[x] + eta*drm->off[x] + drm->on[x+1];

			/* off regime */
			if (x==0) /* lower boundary */
				dnoff[0] = -(1-eta)*r*drm->off[0] - eta*drm->off[0] + drm->on[0];
			else /* off queue */
				dnoff[x] = -(1-eta)*r*drm->off[x] - eta*drm->off[x] + (1-eta)*r*drm->off[x-1];
		}

		/* curtailment mode */
		else
		{
			/* on regime */
			if (x==L-1) /* upper boundary */
				dnon[x] = -drm->on[x] + r*drm->off[x];
			else /* on queue */
				dnon[x] = -drm->on[x] + (1+eta)*drm->on[x+1];

			/* off regime */
			if (x==0) /* lower boundary */
				dnoff[0] = -r*drm->off[0] + drm->on[0];
			else /* off queue */
				dnoff[x] = -r*drm->off[x] - eta*drm->on[x] + r*drm->off[x-1];
		}

		/* check for problem */
		if (dnon[x]+drm->on[x]<0 || dnon[x]+drm->on[x]>1 || dnoff[x]+drm->off[x]<0 || dnoff[x]+drm->off[x]>1)
		{
			errmsg = "drm_update(): model saturated";
			return QNAN;
		}
	}

	/* update states */
	for (x=0; x<L; x++)
	{
		non += drm->on[x] += dnon[x];
		noff += drm->off[x] += dnoff[x];
	}

	return non;
}

/* Calculate the quantity of demand

	Returns:
		The current fractional load.
*/
double drm_quantity(DRMODEL *drm)
{
	double non = 0.0;
	int x;
	for (x=0; x<drm->nbins; x++)
		non += drm->on[x];
	return non;
}

/* Estimate the price of a response

	Arguments:
		DRMODEL *drm		the demand response model
		double Q			the quantity for which the price is desired (between 0 and 1)
		double phi			the natural duty cycle (between 0 and 1)
		double Ed			the elasticity of demand (must be stricly negative)

	Returns:
		The price needed to achieve the change in quantity desired.

	Remarks:
		The demand curve is assumed to be a sigmoid function with
		with the expectation quantity at the mean price and no
		standard deviation.  Increasing load is associated with
		lower prices and decreasing load is associated with higher
		prices.
*/
double drm_price(DRMODEL *drm, 
				 double Q,
				 double phi,
				 double Ed)
{
#define PRICE(Q,D,E) ((log(pow((Q),-(D))-1)-log(D))/-(E))
	const double QNAN = *(double*)nan;
	double Qm, Pm, Pa, nu;

	/* check inputs */
	if (Ed>=0)
	{
		errmsg = "drm_price(): eta must be stricly negative";
		return QNAN;
	}
	if (Q<=0 || Q>=1)
	{
		errmsg = "drm_price(): Q must be exclusively between 0 and 1 ";
		return QNAN;
	}
	if (phi<=0 || phi>=1)
	{
		errmsg = "drm_price(): phi must be exclusively between 0 and 1 ";
		return QNAN;
	}
	nu = 1/(1-(phi-0.5)/0.5);

	/* compute the current load */
	Qm = drm_quantity(drm);

	/* compute the current price */
	Pm = PRICE(Qm,nu,Ed);

	/* compute the objective price */
	Pa = PRICE(Q,nu,Ed);

	return Pa-Pm;
}

/* Estimate the quantity deviation expected at a given price deviation

	Arguments:
		DRMODEL *drm	the demand response model
		double dP		the price deviation to consider
		double Ed		the elasticity of demand (always negative)
		double delta	the band control signal
		double eta		the curtailment control signal
		double phi		the duty cycle control signal

	Returns:
		The expected demand given the conditions given

	Remarks:
		The model is not updated but copies are tested
		until a load fraction is found to with a precision of 0.001
		that has the given price. If no result is found, QNAN is
		returned.
*/
double drm_forecast(DRMODEL *drm, double Px, double Ed, short delta, double eta, double phi)
{
	const double QNAN = *(double*)nan;

	/* make a working copy of the original model */
	DRMODEL *test = drm_copy(drm);
	
	/* get the baseline load */
	double Q0 = drm_update(test,delta,eta,phi);

	/* get the baseline price */
	double P = drm_price(test,Q0,phi,Ed);

	/* get the search bounds */
	double Q1 = Q0;
	double Q2 = Px<P ? 1.0 : 0.0;

	/* next guess */
	double Q = Q0;

	/* begin binary search */
	while (fabs(Q1-Q2)>0.001)
	{
		/* update next guess */
		Q = (Q1+Q2)/2;

		/* calculate the price for this guess */
		P =	drm_price(test,Q,phi,Ed);

		/* update the search bounds based on price */
		if (P<Px)
			Q1 = Q;
		else
			Q2 = Q;
	}

	/* free copy */
	free(test);
	return Q;
}