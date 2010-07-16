/* dr_model.h
 * Copyright (C) 2010 Battelle Memorial Institute
 * All Rights Reserved
 *
 * Author: David Chassin, Pacific Northwest National Laboratory
 *
 */

#ifndef _DRMODEL_H
#define _DRMODEL_H

typedef struct {
	double *on;
	double *off;
	unsigned short nbins;
} DRMODEL;

char *drm_error(void);
DRMODEL *drm_create(unsigned short L, double eta, double phi);
DRMODEL *drm_copy(DRMODEL *from);
int drm_equal(DRMODEL *to, DRMODEL *from);
double drm_update(DRMODEL *drm, short delta, double eta, double phi);
double drm_quantity(DRMODEL *drm);
double drm_price(DRMODEL *drm, double Q, double phi, double Ed);
double drm_forecast(DRMODEL *drm, double Px, double Ed, short delta, double eta, double phi);
double drm_entropy(DRMODEL *drm);

#endif
