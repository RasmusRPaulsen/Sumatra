﻿/**********************************************************************

polygonizer.h

This is Jules Bloomenthal's implicit surface polygonizer from GRAPHICS 
GEMS IV. Bloomenthal's polygonizer is still used and the present code
is simply the original code morphed into C++.

J. Andreas Bærentzen 2003.

**********************************************************************/

#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <vector>
#include <list>
#include <sys/types.h>
#include "BloomenthalPolygonizer.h"

#include <vtkPolyData.h>
#include <vtkPoints.h>
#include <vtkCellArray.h>

using namespace std;

namespace Geometry
{

	namespace
	{
		const int RES =	10; /* # converge iterations    */

		const int L =	0;  /* left direction:	-x, -i */
		const int R =	1;  /* right direction:	+x, +i */
		const int B =	2;  /* bottom direction: -y, -j */
		const int T =	3;  /* top direction:	+y, +j */
		const int N =	4;  /* near direction:	-z, -k */
		const int F =	5;  /* far direction:	+z, +k */
		const int LBN =	0;  /* left bottom near corner  */
		const int LBF =	1;  /* left bottom far corner   */
		const int LTN =	2;  /* left top near corner     */
		const int LTF =	3;  /* left top far corner      */
		const int RBN =	4;  /* right bottom near corner */
		const int RBF =	5;  /* right bottom far corner  */
		const int RTN =	6;  /* right top near corner    */
		const int RTF =	7;  /* right top far corner     */


		/* the LBN corner of cube (i, j, k), corresponds with location
		 * (start.x+(i-.5)*size, start.y+(j-.5)*size, start.z+(k-.5)*size) */

		inline double RAND() 
		{
			return (rand()&32767)/32767.0f;
		}

		const int HASHBIT = 5;

		const int HASHSIZE = (size_t)(1<<(3*HASHBIT));   

		const int MASK = ((1<<HASHBIT)-1);

		inline int HASH( int i, int j,int k) 
		{ 
			return (((((i&MASK)<<HASHBIT)|j&MASK)<<HASHBIT)|k&MASK);
		} 

		inline int BIT(int i, int bit) 
		{ 
			return (i>>bit)&1; 
		}

		// flip the given bit of i 
		inline int FLIP(int i, int bit) 
		{
			return i^1<<bit;
		}

		struct TEST {		   /* test the function for a signed value */
			POINT p;			   /* location of test */
			double value;		   /* function value at p */
			int ok;			   /* if value is of correct sign */
		};

		struct CORNER {		   /* corner of a cube */
			int i, j, k;		   /* (i, j, k) is index within lattice */
			double x, y, z, value;	   /* location and function value */
		};

		struct CUBE {		   /* partitioning cell (cube) */
			int i, j, k;		   /* lattice location of cube */
			CORNER *corners[8];		   /* eight corners */
		};

		struct CENTERELEMENT {	   /* list of cube locations */
			int i, j, k;		   /* cube location */
			CENTERELEMENT(int _i, int _j, int _k): i(_i), j(_j), k(_k) {}
		};
		typedef list<CENTERELEMENT> CENTERLIST;

		struct CORNERELEMENT {	   /* list of corners */
			int i, j, k;		   /* corner id */
			double value;		   /* corner value */
			CORNERELEMENT(int _i, int _j, int _k, double _value):
				i(_i), j(_j), k(_k), value(_value) {}
		};
		typedef list<CORNERELEMENT> CORNERLIST;

		struct EDGEELEMENT {	   /* list of edges */
			int i1, j1, k1, i2, j2, k2;	   /* edge corner ids */
			int vid;			   /* vertex id */
		};

		typedef list<EDGEELEMENT> EDGELIST;
		typedef list<int> INTLIST;
		typedef list<INTLIST> INTLISTS;


		//----------------------------------------------------------------------
		// Implicit surface evaluation functions
		//----------------------------------------------------------------------

		/* converge: from two points of differing sign, converge to zero crossing */

		void converge (POINT* p1, POINT* p2, double v, 
									 vtkImplicitFunction* function, POINT* p)
		{
			int i = 0;
			POINT pos, neg;
			if (v < 0) {
				pos.x = p2->x; pos.y = p2->y; pos.z = p2->z;
				neg.x = p1->x; neg.y = p1->y; neg.z = p1->z;
			}
			else {
				pos.x = p1->x; pos.y = p1->y; pos.z = p1->z;
				neg.x = p2->x; neg.y = p2->y; neg.z = p2->z;
			}
			while (1) {
				p->x = 0.5f*(pos.x + neg.x);
				p->y = 0.5f*(pos.y + neg.y);
				p->z = 0.5f*(pos.z + neg.z);
				if (i++ == RES) return;
				if ((function->FunctionValue(p->x, p->y, p->z)) > 0.0)
					{pos.x = p->x; pos.y = p->y; pos.z = p->z;}
				else {neg.x = p->x; neg.y = p->y; neg.z = p->z;}
			}
		}


		/* vnormal: compute unit length surface normal at point */

		inline void vnormal (vtkImplicitFunction* function, POINT* point, POINT* n, double delta)
		{
			double f = function->FunctionValue(point->x, point->y, point->z);
			n->x = function->FunctionValue(point->x+delta, point->y, point->z)-f;
			n->y = function->FunctionValue(point->x, point->y+delta, point->z)-f;
			n->z = function->FunctionValue(point->x, point->y, point->z+delta)-f;
			f = sqrt(n->x*n->x + n->y*n->y + n->z*n->z);
			if (f != 0.0) {n->x /= f; n->y /= f; n->z /= f;}
		}



		// ----------------------------------------------------------------------

		class EDGETABLE
		{
			vector<EDGELIST> table;		   /* edge and vertex id hash table */
		
		public:
		
			EDGETABLE(): table(2*HASHSIZE) {}

			void setedge (int i1, int j1, int k1, 
										int i2, int j2, int k2, int vid);
			int getedge (int i1, int j1, int k1, 
									 int i2, int j2, int k2);
		

		};

		/* setedge: set vertex id for edge */
		void EDGETABLE::setedge (int i1, int j1, int k1, 
														 int i2, int j2, int k2, int vid)
		{
			unsigned int index;

			if (i1>i2 || (i1==i2 && (j1>j2 || (j1==j2 && k1>k2)))) {
				int t=i1; i1=i2; i2=t; t=j1; j1=j2; j2=t; t=k1; k1=k2; k2=t;
			}
			index = HASH(i1, j1, k1) + HASH(i2, j2, k2);

			EDGEELEMENT new_obj;
			new_obj.i1 = i1; 
			new_obj.j1 = j1; 
			new_obj.k1 = k1;
			new_obj.i2 = i2; 
			new_obj.j2 = j2; 
			new_obj.k2 = k2;
			new_obj.vid = vid;
			table[index].push_front(new_obj);
		}


		/* getedge: return vertex id for edge; return -1 if not set */

		int EDGETABLE::getedge (int i1, int j1, int k1, int i2, int j2, int k2)
		{

			if (i1>i2 || (i1==i2 && (j1>j2 || (j1==j2 && k1>k2)))) {
				int t=i1; i1=i2; i2=t; t=j1; j1=j2; j2=t; t=k1; k1=k2; k2=t;
			};
			int hashval = HASH(i1, j1, k1)+HASH(i2, j2, k2);
			EDGELIST::const_iterator q = table[hashval].begin();
			for (; q != table[hashval].end(); ++q)
				if (q->i1 == i1 && q->j1 == j1 && q->k1 == k1 &&
						q->i2 == i2 && q->j2 == j2 && q->k2 == k2)
					return q->vid;
			return -1;
		}


		// ----------------------------------------------------------------------
		class PROCESS
		{	   /* parameters, function, storage */

			std::vector<NORMAL>* gnormals;  
			std::vector<VERTEX>* gvertices;  
			std::vector<TRIANGLE> *gtriangles;
  
			vtkImplicitFunction* function;	   /* implicit surface function */

			double size, delta;		   /* cube size, normal delta */
			int bounds;			   /* cube range within lattice */
			POINT start;		   /* start point on surface */
			bool use_normals;

			// Global list of corners (keeps track of memory)
			list<CORNER> corner_lst;
			list<CUBE> cubes;		   /* active cubes */
			vector<CENTERLIST> centers;	   /* cube center hash table */
			vector<CORNERLIST> corners;	   /* corner value hash table */
			EDGETABLE edges;

			CORNER *setcorner (int i, int j, int k);

			void testface (int i, int j, int k, CUBE* old, 
										 int face, int c1, int c2, int c3, int c4); 

			TEST find (int sign, double x, double y, double z);
    
			int vertid (CORNER* c1, CORNER* c2);
    
			int dotet (CUBE* cube, int c1, int c2, int c3, int c4);
    
			int docube (CUBE* cube);

			double *VolumeLimits;

			int triangle (int i1, int i2, int i3)
			{
				TRIANGLE t;
				t.v0 = i1;
				t.v1 = i2;
				t.v2 = i3;
				(*gtriangles).push_back(t);
				return 1;
			}

		public:
			PROCESS(vtkImplicitFunction* _function,
							double _size, double _delta, 
							int _bounds,
							vector<VERTEX>& _gvertices,
							vector<NORMAL>& _gnormals,
							vector<TRIANGLE>& _gtriangles,
							double *VolLimits,
							bool _compute_normals=false);
	  

			~PROCESS() {}
		
			void march(int mode, double x, double y, double z);

		};


		/* setcenter: set (i,j,k) entry of table[]
			 return 1 if already set; otherwise, set and return 0 */
		int setcenter(vector<CENTERLIST>& table,
									int i, int j, int k)
		{
			int index = HASH(i,j,k);
			CENTERLIST::const_iterator q = table[index].begin();
			for(; q != table[index].end(); ++q)
				if(q->i == i && q->j==j && q->k == k) return 1;
  
			CENTERELEMENT elem(i,j,k);
			table[index].push_front(elem);
			return 0;
		}

		/* setcorner: return corner with the given lattice location
			 set (and cache) its function value */
		CORNER* PROCESS::setcorner (int i, int j, int k)
		{
			/* for speed, do corner value caching here */
			corner_lst.push_back(CORNER());
			CORNER *c = &corner_lst.back();
			int index = HASH(i, j, k);

			c->i = i; c->x = start.x+((double)i-.5f)*size;
			c->j = j; c->y = start.y+((double)j-.5f)*size;
			c->k = k; c->z = start.z+((double)k-.5f)*size;

			CORNERLIST::const_iterator l = corners[index].begin();
			for (; l != corners[index].end(); ++l)
				if (l->i == i && l->j == j && l->k == k) {
					c->value = l->value;
					return c;
				}

			c->value = function->FunctionValue(c->x, c->y, c->z);
			CORNERELEMENT elem(i,j,k,c->value);
			corners[index].push_front(elem);
			return c;
		}



		/* testface: given cube at lattice (i, j, k), and four corners of face,
		 * if surface crosses face, compute other four corners of adjacent cube
		 * and add new cube to cube stack */

		inline void PROCESS::testface (int i, int j, int k, CUBE* old, 
														int face, int c1, int c2, int c3, int c4) 
		{
			CUBE new_obj;

			static int facebit[6] = {2, 2, 1, 1, 0, 0};
			int n, pos = old->corners[c1]->value > 0.0 ? 1 : 0, bit = facebit[face];

			/* test if no surface crossing, cube out of bounds, or already visited: */
			if ((old->corners[c2]->value > 0) == pos &&
					(old->corners[c3]->value > 0) == pos &&
					(old->corners[c4]->value > 0) == pos) return;
			if (abs(i) > bounds || abs(j) > bounds || abs(k) > bounds) return;
			if (setcenter(centers, i, j, k)) return;

			// \todo The maximum implicit value should be set by the user (Rasmus)
			double maxval = 1e6;

			if (abs(old->corners[c1]->value) > maxval || abs(old->corners[c2]->value) > maxval || abs(old->corners[c3]->value) > maxval || abs(old->corners[c4]->value) > maxval)
				return;

			if (VolumeLimits)
			{
				if (old->corners[c1]->x < VolumeLimits[0] || old->corners[c1]->x > VolumeLimits[1] ||
					old->corners[c1]->y < VolumeLimits[2] || old->corners[c1]->y > VolumeLimits[3] ||
					old->corners[c1]->z < VolumeLimits[4] || old->corners[c1]->z > VolumeLimits[5])
					return;
				if (old->corners[c2]->x < VolumeLimits[0] || old->corners[c2]->x > VolumeLimits[1] ||
					old->corners[c2]->y < VolumeLimits[2] || old->corners[c2]->y > VolumeLimits[3] ||
					old->corners[c2]->z < VolumeLimits[4] || old->corners[c2]->z > VolumeLimits[5])
					return;
				if (old->corners[c3]->x < VolumeLimits[0] || old->corners[c3]->x > VolumeLimits[1] ||
					old->corners[c3]->y < VolumeLimits[2] || old->corners[c3]->y > VolumeLimits[3] ||
					old->corners[c3]->z < VolumeLimits[4] || old->corners[c3]->z > VolumeLimits[5])
					return;
				if (old->corners[c4]->x < VolumeLimits[0] || old->corners[c4]->x > VolumeLimits[1] ||
					old->corners[c4]->y < VolumeLimits[2] || old->corners[c4]->y > VolumeLimits[3] ||
					old->corners[c4]->z < VolumeLimits[4] || old->corners[c4]->z > VolumeLimits[5])
					return;
			}

			/* create new_obj cube: */
			new_obj.i = i;
			new_obj.j = j;
			new_obj.k = k;
			for (n = 0; n < 8; n++) new_obj.corners[n] = 0;
			new_obj.corners[FLIP(c1, bit)] = old->corners[c1];
			new_obj.corners[FLIP(c2, bit)] = old->corners[c2];
			new_obj.corners[FLIP(c3, bit)] = old->corners[c3];
			new_obj.corners[FLIP(c4, bit)] = old->corners[c4];
			for (n = 0; n < 8; n++)
				if (new_obj.corners[n] == 0)
					new_obj.corners[n] = setcorner(i+BIT(n,2), j+BIT(n,1), k+BIT(n,0));

			// Add new cube to top of stack
			cubes.push_front(new_obj);
		}

		/* find: search for point with value of given sign (0: neg, 1: pos) */

		TEST PROCESS::find (int sign, double x, double y, double z)
		{
			int i;
			TEST test;
			double range = size;
			test.ok = 1;
			for (i = 0; i < 10000; i++) {
				test.p.x = x+range*(RAND()-0.5f);
				test.p.y = y+range*(RAND()-0.5f);
				test.p.z = z+range*(RAND()-0.5f);
				test.value = function->FunctionValue(test.p.x, test.p.y, test.p.z);
				if (sign == (test.value > 0.0)) return test;
				range = range*1.0005f; /* slowly expand search outwards */
			}
			test.ok = 0;
			return test;
		}


		/* vertid: return index for vertex on edge:
		 * c1->value and c2->value are presumed of different sign
		 * return saved index if any; else compute vertex and save */

		int PROCESS::vertid (CORNER* c1, CORNER* c2)
		{
			VERTEX v;
			NORMAL n;
			POINT a, b;
			int vid = edges.getedge(c1->i, c1->j, c1->k, c2->i, c2->j, c2->k);
			if (vid != -1) return vid;			     /* previously computed */
			a.x = c1->x; a.y = c1->y; a.z = c1->z;
			b.x = c2->x; b.y = c2->y; b.z = c2->z;
			converge(&a, &b, c1->value, function, &v); /* position */
			(*gvertices).push_back(v);			   /* save vertex */
			if(use_normals)
				{
					vnormal(function, &v, &n, delta);			   /* normal */
					(*gnormals).push_back(n);			   /* save vertex */
				}
			vid = gvertices->size()-1;
			edges.setedge(c1->i, c1->j, c1->k, c2->i, c2->j, c2->k, vid);
			return vid;
		}




		/**** Tetrahedral Polygonization ****/


		/* dotet: triangulate the tetrahedron
		 * b, c, d should appear clockwise when viewed from a
		 * return 0 if client aborts, 1 otherwise */

		int PROCESS::dotet (CUBE* cube, int c1, int c2, int c3, int c4)
		{
			CORNER *a = cube->corners[c1];
			CORNER *b = cube->corners[c2];
			CORNER *c = cube->corners[c3];
			CORNER *d = cube->corners[c4];
			int index = 0, apos, bpos, cpos, dpos, e1, e2, e3, e4, e5, e6;
			if ((apos = (a->value > 0.0))) index += 8;
			if ((bpos = (b->value > 0.0))) index += 4;
			if ((cpos = (c->value > 0.0))) index += 2;
			if ((dpos = (d->value > 0.0))) index += 1;
			/* index is now 4-bit number representing one of the 16 possible cases */
			if (apos != bpos) e1 = vertid(a, b);
			if (apos != cpos) e2 = vertid(a, c);
			if (apos != dpos) e3 = vertid(a, d);
			if (bpos != cpos) e4 = vertid(b, c);
			if (bpos != dpos) e5 = vertid(b, d);
			if (cpos != dpos) e6 = vertid(c, d);
			/* 14 productive tetrahedral cases (0000 and 1111 do not yield polygons */
			switch (index) {
			case 1:	 return triangle(e5, e6, e3);
			case 2:	 return triangle(e2, e6, e4);
			case 3:	 return triangle(e3, e5, e4) &&
								 triangle(e3, e4, e2);
			case 4:	 return triangle(e1, e4, e5);
			case 5:	 return triangle(e3, e1, e4) &&
								 triangle(e3, e4, e6);
			case 6:	 return triangle(e1, e2, e6) &&
								 triangle(e1, e6, e5);
			case 7:	 return triangle(e1, e2, e3);
			case 8:	 return triangle(e1, e3, e2);
			case 9:	 return triangle(e1, e5, e6) &&
								 triangle(e1, e6, e2);
			case 10: return triangle(e1, e3, e6) &&
								 triangle(e1, e6, e4);
			case 11: return triangle(e1, e5, e4);
			case 12: return triangle(e3, e2, e4) &&
								 triangle(e3, e4, e5);
			case 13: return triangle(e6, e2, e4);
			case 14: return triangle(e5, e3, e6);
			}
			return 1;
		}


		/**** Cubical Polygonization (optional) ****/


		const int LB =	0;  /* left bottom edge	*/
		const int LT =	1;  /* left top edge	*/
		const int LN =	2;  /* left near edge	*/
		const int LF =	3;  /* left far edge	*/
		const int RB =	4;  /* right bottom edge */
		const int RT =	5;  /* right top edge	*/
		const int RN =	6;  /* right near edge	*/
		const int RF =	7;  /* right far edge	*/
		const int BN =	8;  /* bottom near edge	*/
		const int BF =	9;  /* bottom far edge	*/
		const int TN =	10; /* top near edge	*/
		const int TF =	11; /* top far edge	*/


		class CUBETABLE
		{
			vector<INTLISTS> ctable;
			int nextcwedge (int edge, int face);
			int otherface (int edge, int face);

		public:

			CUBETABLE();
		
			const INTLISTS& get_lists(int i) const
			{
				return ctable[i];
			}
		};

		/*			edge: LB, LT, LN, LF, RB, RT, RN, RF, BN, BF, TN, TF */
		const int corner1[12]	   = {LBN,LTN,LBN,LBF,RBN,RTN,RBN,RBF,LBN,LBF,LTN,LTF};
		const int corner2[12]	   = {LBF,LTF,LTN,LTF,RBF,RTF,RTN,RTF,RBN,RBF,RTN,RTF};
		const int leftface[12]	   = {B,  L,  L,  F,  R,  T,  N,  R,  N,  B,  T,  F};
		/* face on left when going corner1 to corner2 */
		const int rightface[12]   = {L,  T,  N,  L,  B,  R,  R,  F,  B,  F,  N,  T};
		/* face on right when going corner1 to corner2 */


		/* nextcwedge: return next clockwise edge from given edge around given face */

		inline int CUBETABLE::nextcwedge (int edge, int face)
		{
			switch (edge) {
			case LB: return (face == L)? LF : BN;
			case LT: return (face == L)? LN : TF;
			case LN: return (face == L)? LB : TN;
			case LF: return (face == L)? LT : BF;
			case RB: return (face == R)? RN : BF;
			case RT: return (face == R)? RF : TN;
			case RN: return (face == R)? RT : BN;
			case RF: return (face == R)? RB : TF;
			case BN: return (face == B)? RB : LN;
			case BF: return (face == B)? LB : RF;
			case TN: return (face == T)? LT : RN;
			case TF: return (face == T)? RT : LF;
			}
			return -1;
		}

		/* otherface: return face adjoining edge that is not the given face */

		inline int CUBETABLE::otherface (int edge, int face)
		{
			int other = leftface[edge];
			return face == other? rightface[edge] : other;
		}


		CUBETABLE::CUBETABLE(): ctable(256)
		{
			int i, e, c, done[12], pos[8];
			for (i = 0; i < 256; i++) 
				{
					for (e = 0; e < 12; e++) 
						done[e] = 0;
					for (c = 0; c < 8; c++) 
						pos[c] = BIT(i, c);
					for (e = 0; e < 12; e++)
						if (!done[e] && (pos[corner1[e]] != pos[corner2[e]])) 
							{
								INTLIST ints;
								int start = e, edge = e;
							
								/* get face that is to right of edge from pos to neg corner: */
								int face = pos[corner1[e]]? rightface[e] : leftface[e];
								while (1) 
									{
										edge = nextcwedge(edge, face);
										done[edge] = 1;
										if (pos[corner1[edge]] != pos[corner2[edge]]) 
											{
												ints.push_front(edge);
												if (edge == start) 
													break;
												face = otherface(edge, face);
											}
									}
								ctable[i].push_front(ints);
							}
				}
		}

		inline const INTLISTS& get_cubetable_entry(int i) 
		{
			static CUBETABLE c;
			return c.get_lists(i);
		}


		/* docube: triangulate the cube directly, without decomposition */

		int PROCESS::docube (CUBE* cube)
		{
			int index = 0;
			for (int i = 0; i < 8; i++) 
				if (cube->corners[i]->value > 0.0) 
					index += (1<<i);

			INTLISTS intlists = get_cubetable_entry(index);
			INTLISTS::const_iterator polys = intlists.begin();
			for (; polys != intlists.end(); ++polys) 
				{
					INTLIST::const_iterator edges = polys->begin();
					int a = -1, b = -1, count = 0;
					for (; edges != polys->end(); ++edges) 
						{
							CORNER *c1 = cube->corners[corner1[(*edges)]];
							CORNER *c2 = cube->corners[corner2[(*edges)]];
							int c = vertid(c1, c2);
							if (++count > 2 && ! triangle(a, b, c)) 
								return 0;
							if (count < 3) 
								a = b;
							b = c;
						}
				}
			return 1;
		}







		/**** An Implicit Surface BloomenthalPolygonizer ****/


		/* polygonize: polygonize the implicit surface function
		 *   arguments are:
		 *	 ImplicitFunction
		 *	     the implicit surface function
		 *	     return negative for inside, positive for outside
		 *	 double size
		 *	     width of the partitioning cube
		 *   double delta 
		 *       a small step - used for gradient computation
		 *	 int bounds
		 *	     max. range of cubes (+/- on the three axes) from first cube
		 *   _gvertices, _gnormals, _gtriangles
		 *       the data structures into which information is put.
		 */

		PROCESS::PROCESS(vtkImplicitFunction* _function,
										 double _size, double _delta, 
										 int _bounds, 
										 vector<VERTEX>& _gvertices,
										 vector<NORMAL>& _gnormals,
										 vector<TRIANGLE>& _gtriangles,
										 double *VolLimits,
										 bool _use_normals):
			function(_function), size(_size), delta(_delta), bounds(_bounds),
			centers(HASHSIZE), corners(HASHSIZE), 
			gvertices(&_gvertices),
			gnormals(&_gnormals),
			gtriangles(&_gtriangles),
			VolumeLimits(VolLimits),
			use_normals(_use_normals)
		{}

		void PROCESS::march(int mode, double x, double y, double z)
		{
			int noabort;
			TEST in, out;
  
			/* find point on surface, beginning search at (x, y, z): */
			srand(1);
			in = find(1, x, y, z);
			out = find(0, x, y, z);
			if (!in.ok || !out.ok) 
				throw(string("can't find starting point"));
  
			converge(&in.p, &out.p, in.value, function, &start);
  
			/* push initial cube on stack: */
			CUBE cube;
			cube.i = cube.j = cube.k = 0;
			cubes.push_front(cube);

			/* set corners of initial cube: */
			for (int n = 0; n < 8; n++)
				cubes.front().corners[n] = setcorner(BIT(n,2), BIT(n,1), BIT(n,0));
    
			setcenter(centers, 0, 0, 0);
    
			while (cubes.size() != 0) 
				{
					/* process active cubes till none left */
					CUBE c = cubes.front();
      
					noabort = mode == TET?
						/* either decompose into tetrahedra and polygonize: */
						dotet(&c, LBN, LTN, RBN, LBF) &&
						dotet(&c, RTN, LTN, LBF, RBN) &&
						dotet(&c, RTN, LTN, LTF, LBF) &&
						dotet(&c, RTN, RBN, LBF, RBF) &&
						dotet(&c, RTN, LBF, LTF, RBF) &&
						dotet(&c, RTN, LTF, RTF, RBF)
						:
						/* or polygonize the cube directly: */
						docube(&c);
					if (! noabort) throw string("aborted");
      
					/* pop current cube from stack */
					cubes.pop_front();
      
					/* test six face directions, maybe add to stack: */
					testface(c.i-1, c.j, c.k, &c, L, LBN, LBF, LTN, LTF);
					testface(c.i+1, c.j, c.k, &c, R, RBN, RBF, RTN, RTF);
					testface(c.i, c.j-1, c.k, &c, B, LBN, LBF, RBN, RBF);
					testface(c.i, c.j+1, c.k, &c, T, LTN, LTF, RTN, RTF);
					testface(c.i, c.j, c.k-1, &c, N, LBN, LTN, RBN, RTN);
					testface(c.i, c.j, c.k+1, &c, F, LBF, LTF, RBF, RTF);
				}
		}
	}

	void BloomenthalPolygonizer::march(double x, double y, double z)
		{
			gvertices.clear();
			gnormals.clear();
			gtriangles.clear();
			PROCESS p(func, size, size/(double)(RES*RES), bounds, 
								gvertices, gnormals, gtriangles, VolumeLimits, use_normals);
			p.march(use_tetra?TET:NOTET,x,y,z);
		}

	//! Return surface as a VTK polydata
	void BloomenthalPolygonizer::ExportToPolyData(vtkPolyData *pd)
	{
		int i;

		vtkPoints *pts = vtkPoints::New();
		vtkCellArray *newPolys = vtkCellArray::New();

		int nverts = no_vertices();

		for (i = 0; i < nverts; i++)
		{
			double p[3];
			p[0] = gvertices[i].x;
			p[1] = gvertices[i].y;
			p[2] = gvertices[i].z;

			vtkIdType id = pts->InsertNextPoint(p);
		}
	
		int ntriangles = no_triangles();

		for (int f = 0; f < ntriangles; f++)
		{
			newPolys->InsertNextCell(3);

			vtkIdType id1 = gtriangles[f].v0;
			newPolys->InsertCellPoint(id1);
			vtkIdType id2 = gtriangles[f].v1;
			newPolys->InsertCellPoint(id2);
			vtkIdType id3 = gtriangles[f].v2;
			newPolys->InsertCellPoint(id3);
		}

		pd->SetPoints(pts);
		pd->SetPolys(newPolys);
		newPolys->Delete();
		pts->Delete();
	}
}	

