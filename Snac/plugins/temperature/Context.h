/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
**
** Copyright (C), 2003, 
**	Steve Quenette, 110 Victoria Street, Melbourne, Victoria, 3053, Australia.
**	Californian Institute of Technology, 1200 East California Boulevard, Pasadena, California, 91125, USA.
**	University of Texas, 1 University Station, Austin, Texas, 78712, USA.
**
** Authors:
**	Stevan M. Quenette, Senior Software Engineer, VPAC. (steve@vpac.org)
**	Stevan M. Quenette, Visitor in Geophysics, Caltech.
**	Luc Lavier, Research Scientist, The University of Texas. (luc@utig.ug.utexas.edu)
**	Luc Lavier, Research Scientist, Caltech.
**
** This program is free software; you can redistribute it and/or modify it
** under the terms of the GNU General Public License as published by the
** Free Software Foundation; either version 2, or (at your option) any
** later version.
** 
** This program is distributed in the hope that it will be useful,
** but WITHOUT ANY WARRANTY; without even the implied warranty of
** MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
** GNU General Public License for more details.
** 
** You should have received a copy of the GNU General Public License
** along with this program; if not, write to the Free Software
** Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
**
*/
/** \file
** Role:
**
** Assumptions:
**
** Comments:
**
** $Id: Context.h 3265 2006-11-15 21:43:39Z EunseoChoi $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#ifndef __SnacTemperature_Context_h__
#define __SnacTemperature_Context_h__
	
	/* Context Information */
	struct _SnacTemperature_Context {
		/* For _SnacTemperature_Top2BottomSweep condition function */
		double				topTemp;
		double				bottomTemp;
		
	  //for SnacTemperature_erf
	  double v_stretch;
	  double startX;
	  double endX;
	  double                    crustal_thickness;
	  double                    crustal_thermal_gradient;
	  double                    T1; //first layer bottom temperature or erf top temperature at 5km depth  

	  double      minY;

	  //For SnacTemperature_Transform_Fault
	  double Segment1_x;
	  double Segment1_z_min;
	  double Segment1_z_max;
	  
	  double Segment2_x;
	  double Segment2_z_min;
	  double Segment2_z_max;
	  
		/* ICs and BCs */
		CompositeVC*			temperatureBCs;
		
		/* Dumping */
		FILE*				temperatureOut;
		FILE*				temperatureCheckpoint;
	};
	
	/* Print the contents of the context extension */
	void SnacTemperature_Context_Print( void* _context, Stream* stream );
	
#endif /* __SnacTemperature_Context_h__ */
