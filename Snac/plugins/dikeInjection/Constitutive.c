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
** $Id: Constitutive.c 3274 2007-03-27 20:25:29Z EunseoChoi $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StGermain/FD/FD.h>
#include "Snac/Snac.h"
#include "Snac/ViscoPlastic/ViscoPlastic.h"
#include "types.h"
#include "Context.h"
#include "Constitutive.h"
#include "Register.h"
#include <math.h>
#include <string.h>
#include <assert.h>

#ifndef PI
#ifndef M_PIl
#ifndef M_PI
#define PI 3.14159265358979323846
#else
#define PI M_PI
#endif
#else
#define PI M_PIl
#endif
#endif
//#define DEBUG

//int maxxK = 9;  //meshSizeK - 2 ; for M-factor vary from 0 to 1

void SnacDikeInjection_Constitutive( void* _context, Element_LocalIndex element_lI ) {

      
	Snac_Context*		        context = (Snac_Context*)_context;
	SnacDikeInjection_Context*	contextExt = ExtensionManager_Get(
					       		context->extensionMgr,
	      						context,
							SnacDikeInjection_ContextHandle );
	Snac_Element			*element = Snac_Element_At( context, element_lI );
	SnacViscoPlastic_Element* viscoplasticElement = ExtensionManager_Get( context->mesh->elementExtensionMgr, element, SnacViscoPlastic_ElementHandle );
	const Snac_Material		*material = &context->materialProperty[element->material_I];
	
	/* make local copies. */
	double startX = contextExt->startX;
	double endX = contextExt->endX;
	double startZ = contextExt->startZ;
	double endZ = contextExt->endZ;
	double dX = endX-startX;
	double dZ = endZ-startZ;
	double Mb = contextExt->Mb;
        double Me = contextExt->Me; 
	//double maxK = contextExt->maxK; //new added;   why this not working        
	//	fprintf (stderr, "maxK1 = %f\n", maxK);
        double elem_dX = 0.0;
	double epsilon_xx = 0.0;
	Tetrahedra_Index tetra_I;//build/include/Snac/types.h:typedef unsigned intTetrahedra_Index;

	/* Some convenience stuffs. */
        //MeshLayout*			meshLayout = (MeshLayout*)context->meshLayout;
	//HexaMD*				decomp = (HexaMD*)meshLayout->decomp;
	IJK				ijk;
	

	//	const Index                     kCount = decomp->nodeGlobal3DCounts[2];
	Mesh*                           mesh = context->mesh;
        MeshLayout*                     layout = (MeshLayout*)mesh->layout;
        HexaMD*                         decomp = (HexaMD*)layout->decomp;

	Element_GlobalIndex             global_K_range = decomp->elementGlobal3DCounts[2];
	
        Element_GlobalIndex		element_gI = _MeshDecomp_Element_LocalToGlobal1D( decomp, element_lI );
        RegularMeshUtils_Element_1DTo3D( decomp, element_gI, &ijk[0], &ijk[1], &ijk[2] );

        //iiiii = N;   why this doesn't work??????????
	//why this doesn't work? 12345678910 012345678910 instead of 1 to n        int iiiii = 0;   //added later for M-factor variation
       
	/* for ( int iii = 0; iii<7; iii++) {
	fprintf(stderr," (Snac_Element_NodeCoord( context, element_lI, %d))[0] = %e\n", iii,  Snac_Element_NodeCoord( context, element_lI, iii)[0]); 
	}*/
	//epsilon_xx = (contextExt->injectionRate*context->dt)/elem_dX;

	//fprintf(stderr,"context->dt=%e epsilon_xx=%e\n, injectionRate=%e",context->dt,epsilon_xx, contextExt->injectionRate);
	//	epsilon_xx = (contextExt->injectionRate*context->dt)/elem_dX; 
	
	for( tetra_I = 0; tetra_I < Tetrahedra_Count; tetra_I++ ) {
		Coord baryCenter;
		double distance = 0.0;
		double numer = 0.0;
		double denom = 1.0;
		Node_LocalIndex node_lI;
		unsigned int dim;

	
	

		//                             iiiii = &kk;  	


		//fprintf(stderr, "Tetrahedra_Count=%d,  tetra_I=%d\n",Tetrahedra_Count,tetra_I);

		/*
		  First decide whther this tet is a part of the dike.
		*/

		/* compute barycenter. */
		for( dim=0; dim < 3; dim++ ){
		   baryCenter[dim] = 0.0;
		}

               for(node_lI=0; node_lI<4; node_lI++) {
	     Coord* tetNodeCoord = Snac_Element_NodeCoord( context, element_lI, TetraToNode[tetra_I][node_lI] );
	     //fprintf(stderr, "tetra_I=%d\n", tetra_I);
	     //fprintf(stderr, "node_lI=%d\n", node_lI);
	     //fprintf(stderr, "tetNodeCoord[%d][%d]=%e\n", tetra_I, node_lI, Snac_Element_NodeCoord( context, element_lI, TetraToNode[tetra_I][node_lI]));
	         for( dim=0; dim < 3; dim++ ){
                   baryCenter[dim] += 0.25 * (*tetNodeCoord)[dim];
	           //fprintf(stderr, "(*tetNodeCoord)[%d]=%e\n", dim, (*tetNodeCoord)[dim]);
	         }
               }
	     //fprintf(stderr,"baryCenter[0]=%e\n,baryCenter[1]=%e\n,baryCenter[2]=%e\n",baryCenter[0],baryCenter[1],baryCenter[2]);
        	
/* The following is the general formula for distance from a li（）(diagonal of dx dz)o a point. */

/*		numer = fabs( dX*(startZ-baryCenter[2])-(startX-baryCenter[0])*dZ );
		denom = sqrt( dX*dX + dZ*dZ );
		assert( denom > 0.0 );
		distance = numer/denom;
*/	
	//fprintf(stderr, "distance = %e\n", distance);	
          	//	fprintf(stderr,"distance=%e, dikeDepth=%e",distance,contextExt->dikeDepth);
	       
/*the above way for restricting dike in x-axis can be replaced by simple method be low*/		
	       distance = fabs( baryCenter[0] - 0.5 * (startX + endX) );

                 /* 
		   If part of the dike, adjust stresses. 

		   Note that although parameters can define a ridge with an arbitrary orientation,
		   the following stress mods assume ridges are parallel with z-axis.
		   One can implement tensor rotation when necessary.
		 */
	      

             if( (distance < 0.5*contextExt->dikeWidth) && (baryCenter[1] >= -contextExt->dikeDepth) ){

	elem_dX = 0.25*( 
      	(Snac_Element_NodeCoord( context, element_lI, 1)[0]-Snac_Element_NodeCoord( context, element_lI, 0)[0]) + 
	(Snac_Element_NodeCoord( context, element_lI, 2)[0]-Snac_Element_NodeCoord( context, element_lI, 3)[0]) + 
	(Snac_Element_NodeCoord( context, element_lI, 5)[0]-Snac_Element_NodeCoord( context, element_lI, 4)[0]) + 
        (Snac_Element_NodeCoord( context, element_lI, 6)[0]-Snac_Element_NodeCoord( context, element_lI, 7)[0]) 
	 );
	
	//fprintf(stderr,"elem_dX=%e \n",elem_dX);

	epsilon_xx = (contextExt->injectionRate*context->dt)/elem_dX/2; 

	//fprintf(stderr,"epsilon_xx=%e \n",epsilon_xx);
	         StressTensor*		stress = &element->tetra[tetra_I].stress;
			//fprintf(stderr, "baryCenter[1]=%e\n", baryCenter[1]);
            		//fprintf(stderr, "distance = %e\n", distance);
                        //fprintf(stderr, "element_lI=%d\n", element_lI);				
// fprintf(stderr,"el=%d (%d %d %d) tet=%d (%e %e %e) distance=%e startX=%e width=%e\n",	element_lI,ijk[0],ijk[1],ijk[2],tetra_I,baryCenter[0], baryCenter[1], baryCenter[2],distance,startX,contextExt->dikeWidth);
		 //fprintf(stderr, "origin(*stress)[0][0]=%e\n", (*stress)[0][0]);
/*P197 Geodynamics*/
		 //fprintf(stderr, "ijk[2]=%d\n", ijk[2]);                 
                 int TT = ijk[2];
		 // fprintf(stderr, "maxK = %f\n", maxK);
	      

		 /*  why this way , maxx K will be hardwire to 0???????????????????????????????????????????????
              (*stress)[0][0] -= (material->lambda + 2.0f * material->mu) * epsilon_xx * ijk[2] / maxxK;
              (*stress)[1][1] -= material->lambda * epsilon_xx * ijk[2] / maxxK;
              (*stress)[2][2] -= material->lambda * epsilon_xx * ijk[2] / maxxK;
		 */
		 //the if global_K_range==1 is for 1d pseudo test
		 if (global_K_range==1){
		 (*stress)[0][0] -= (material->lambda + 2.0f * material->mu) * epsilon_xx;
		 (*stress)[1][1] -= material->lambda * epsilon_xx;
		 (*stress)[2][2] -= material->lambda * epsilon_xx;
		 
                 }else{
		   (*stress)[0][0] -= (material->lambda + 2.0f * material->mu) * epsilon_xx * ((ijk[2]+1) / ( global_K_range)*(Me-Mb)+Mb);
		 (*stress)[1][1] -= material->lambda * epsilon_xx * ((ijk[2]+1) / ( global_K_range)*(Me-Mb)+Mb);
		 (*stress)[2][2] -= material->lambda * epsilon_xx * ((ijk[2]+1) / ( global_K_range)*(Me-Mb)+Mb);
		 }		 
//fprintf(stderr, " global_K_range=%d\n ijk[2]=%d\n",  global_K_range, ijk[2]);
		 
			//fprintf(stderr, "dike(*stress)[0][0]=%e\n", (*stress)[0][0]);			
			/* Also assuming viscoplastic rheology is used. */
		 //		       	viscoplasticElement->plasticStrain[tetra_I] = 0.0;
			//just add the above line code can solve the plastic strain too much at dike problem!!
			//However, I took 4hours playing around in viscoplastic_BI
			//try to add everything here at there, which is so hard
			//adding a variable i.e. endZ itsel is already hard, let along finding the diking element and 
			//get zero plasticStrain
			//Point is
			//1)Before doing thing, really think carefully, sometimes life will be much easier
			// when you use a easier method
			//2) Pointer is great, it can let things easier.
			//3)you have to know exactly what is the global variable, what does it contains, which 
			//pointers can you use everywhere and how to use it. 
	     }

             if( (distance < contextExt->dikeWidth) && (baryCenter[1] >= -contextExt->dikeDepth) ){

 	viscoplasticElement->plasticStrain[tetra_I] = 0.0;


	     }
	}
}
