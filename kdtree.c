#define _GNU_SOURCE

#include <string.h>
#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>
#include <stddef.h>
#include <math.h>
#include <ctype.h>
#include <float.h>

#include "kdtree.h"

/**
 * A set of geographic locations given by latitude and longitude
 * implemented using a k-d tree, where k = 2.  For all functions,
 * points are considered different based on the values of the
 * coordinates that are passed in, even if the coordinates represent
 * the same physical point, so, for example, longitude -180 is
 * considered different than longitude 180.
 */

typedef struct node node;

int index_in_sub(location *list, int list_length, bool x, double value);
location *balance_locations(location *list, location *other_half, int start, int end, int other_half_length, int length, int gen, location *out, bool x_sub_arr);
// location *balance_locations(location *list, int length, int gen, location *out);
location *loc_add(location *list, location *addendum);
void print_list(location *list, int length, bool for_x);
location *sub_arr(location *list, int start, int end);
int median_index(int length);
void mergeSort(int n, const location a[], location out[], bool x);
void mergex(int n1, const location a1[], int n2, const location a2[], location out[]);
void mergey(int n1, const location a1[], int n2, const location a2[], location out[]);

bool loc_equals(location *l1, location *l2);
int kdtree_add_aux(node *head, const location *p, int tree_depth, int count);

bool kdtree_contains_aux(node *head, const location *p);

void produce_regions(node *head, double BLx, double BLy, double TRx, double TRy, int gen, bool is_head);
void print_region(double A, double B, double BLx, double BLy, double TRx, double TRy, int gen);
node *starting_point(node *head, location *p, int count, int depth);
node *find_neighbor(node *head, location *p, location *curr_neighbor, double *curr_min_dist, bool is_head, node *out);

void traverse(node *head, void (*f)(const location *, void *), void *arg);
void print_root(const location *root, void *arg);
void print_kdtree(kdtree *t);

void traverse_free(node *head);

struct node
{
	location *root;
	node *left;
	node *right;
	double *region;
};


struct _kdtree
{
	node *head;
	int depth;
};



/**
 * Creates a balanced k-d tree containing copies of the points in the
 * given array of locations.  If the array is NULL and n is 0 then the
 * returned tree is empty.  If the tree could not be created then the
 * returned value is NULL.
 *
 * @param pts an array of unique valid locations, non-NULL if n != 0
 * @param n the number of points in that array, or 0 if pts is NULL
 * @return a pointer to the newly created collection of points, or NULL
 * if it could not be created
 */
kdtree *kdtree_create(const location *pts, int n)
{
	kdtree *tree = calloc(1, sizeof(kdtree));
	if (tree == NULL)
	{
		return NULL;
	}

	if (pts == NULL || n == 0)
	{
		tree->head = NULL;
		tree->depth = 0;
	}
	else
	{
		location *pts_copy = malloc(n * sizeof(location));
		for (int i = 0; i < n; i++)
		{
			pts_copy[i].x = pts[i].x;
			pts_copy[i].y = pts[i].y;
		}
		location *out = malloc (n * sizeof(location));
		for (int i = 0; i < n; i++)
		{
			out[i].x = out[i].y = DBL_MAX;
		}
		location *pts_balanced = balance_locations(pts_copy, NULL, 0, n, n, n, 0, out, true);

		free(pts_copy);
		// print_list(pts_balanced, n, true);

		for (int i = 0; i < n; i++)
		{
			if (!kdtree_add(tree, pts_balanced + i))
			{						
				free(out);
				return NULL;
			}
		}
		free(out);
	}

	return tree;
}

void print_list(location *list, int length, bool for_x)
{
	if (for_x)
	{
		printf("X: ");
	}
	else
	{
		printf("Y: ");
	}
	printf("[");
	for (int i = 0; i < length; i++)
	{
		printf("(%lf, %lf)", list[i].x, list[i].y);
		if (i != length - 1)
		{
			printf(", ");
		}
	}
	printf("]\n");
}

int median_index(int length)
{
	return (int) (ceil ((double) length / 2.0)) - 1;
}

/* merge sorted arrays a1 and a2, putting result in out */
void mergex(int n1, const location a1[], int n2, const location a2[], location out[])
{
    int i1;
    int i2;
    int iout;

    i1 = i2 = iout = 0;

    while(i1 < n1 || i2 < n2) {
        if(i2 >= n2 || ((i1 < n1) && (a1[i1].x < a2[i2].x))) {
            /* a1[i1] exists and is smaller */
            out[iout++] = a1[i1++];
        } 
        else if (i2 >= n2 || ((i1 < n1) && (a1[i1].x == a2[i2].x)))
        {
        	if (a1[i1].y < a2[i2].y)
        	{
             	out[iout++] = a1[i1++];
        	}
        	else
        	{
        		out[iout++] = a2[i2++];
        	}
        }


         else {
            /* a1[i1] doesn't exist, or is bigger than a2[i2] */
            out[iout++] = a2[i2++];
        }
    }
}


void mergey(int n1, const location a1[], int n2, const location a2[], location out[])
{
    int i1;
    int i2;
    int iout;

    i1 = i2 = iout = 0;

    while(i1 < n1 || i2 < n2) {
        if(i2 >= n2 || ((i1 < n1) && (a1[i1].y < a2[i2].y))) {
            /* a1[i1] exists and is smaller */
            out[iout++] = a1[i1++];
        } 
        else if (i2 >= n2 || ((i1 < n1) && (a1[i1].y == a2[i2].y)))
        {
        	if (a1[i1].x < a2[i2].x)
        	{
             	out[iout++] = a1[i1++];
        	}
        	else
        	{
        		out[iout++] = a2[i2++];
        	}
        }

         else {
            /* a1[i1] doesn't exist, or is bigger than a2[i2] */
            out[iout++] = a2[i2++];
        }
    }
}

/* sort a, putting result in out */
/* we call this mergeSort to avoid conflict with mergesort in libc */
void mergeSort(int n, const location a[], location out[], bool x)
{
    location *a1;
    location *a2;

    if(n < 2) {
        /* 0 or 1 elements is already sorted */
        memcpy(out, a, sizeof(location) * n);
    } else {
        /* sort into temp arrays */
        a1 = malloc(sizeof(location) * (n/2));
        a2 = malloc(sizeof(location) * (n - n/2));

        if (x)
        {
	        mergeSort(n/2, a, a1, true);
	        mergeSort(n - n/2, a + n/2, a2, true);
	    }
	    else
	    {
	    	mergeSort(n/2, a, a1, false);
	        mergeSort(n - n/2, a + n/2, a2, false);
	    }

        /* merge results */
        if (x)
        {
        	mergex(n/2, a1, n - n/2, a2, out);
        }
        else
        {
        	mergey(n/2, a1, n - n/2, a2, out);
        }

        /* free the temp arrays */
        free(a1);
        free(a2);
    }
}

location *sub_arr(location *list, int start, int end)
{
	location *sub_arr = malloc((end - start)*sizeof(location));
	int count = 0;
	for (int i = start; i < end; i++)
	{
		sub_arr[count] = list[i];
		count++;
	}
	return sub_arr;
}

location *loc_add(location *list, location *addendum)
{
	int i = 0;
	while (list[i].x != DBL_MAX)
	{
		i++;
	}
	list[i] = *addendum;
	return list;
}



int index_in_sub(location *list, int list_length, bool x, double value)
{
	for (int i = 0; i < list_length; i++)
	{
		if (x)
		{
			if (list[i].x == value)
			{
				return i;
			}
		}
		else
		{
			if (list[i].y == value)
			{
				return i;
			}
		}
	}
	return -1;
}

location *balance_locations(location *list, location *other_half, int start, int end, int other_half_length, int length, int gen, location *out, bool x_sub_arr)
{
	// Building the balanced tree
	if (length > 0)
	{
		// creating by_x:
		location *by_x = malloc(length * sizeof(location));

		if (gen >= 0)
		{
			location *list_copy = malloc(length * sizeof(location));
			int count = 0;
			for (int i = start; i < end; i++)
			{
				list_copy[count].x = list[i].x;
				list_copy[count].y = list[i].y;
				count++;
			}
			mergeSort(length, list_copy, by_x, true);
			free(list_copy);
		}
		else
		{
			if (x_sub_arr)
			{
				int count = 0;
				for (int i = start; i < end; i++)
				{
					by_x[count].x = list[i].x;
					by_x[count].y = list[i].y;
					count++;
				}
			}
			else // sorting a Y sub arr by X
			{
				int count = 0; // number of points added in new by_x
				for (int i = 0; i < other_half_length; i++)
				{
					if ((other_half[i].y >= list[start].y && other_half[i].y <= list[end - 1].y) ||
						(other_half[i].y == list[start].y && other_half[i].x >= list[start].x && other_half[i].y <= list[end - 1].y))
					{
						by_x[count].x = other_half[i].x;
						by_x[count].y = other_half[i].y;
						count++;
					}
				}
			}
		}
		// print_list(by_x, length, true);
		// at this point you have sorted by_x

		// print_list(by_x, length, true);
		//

		// creating by_y: 
		location *by_y = malloc(length * sizeof(location));

		if (gen >= 0)
		{
			location *list_copy = malloc(length * sizeof(location));
			int count = 0;
			for (int i = start; i < end; i++)
			{
				list_copy[count].x = list[i].x;
				list_copy[count].y = list[i].y;
				count++;
			}
			mergeSort(length, list_copy, by_y, false);
			free(list_copy);
		}
		else
		{
			if (!x_sub_arr)
			{
				int count = 0;
				for (int i = start; i < end; i++)
				{
					by_y[count].x = list[i].x;
					by_y[count].y = list[i].y;
					count++;
				}			
			}
			else // sorting an X sub arr by Y
			{
				int count = 0; // number of points added in new by_y
				for (int i = 0; i < other_half_length; i++)
				{
					if ((other_half[i].x > list[start].x && other_half[i].x <= list[end - 1].x) ||
						(other_half[i].x == list[start].x && other_half[i].y >= list[start].y && other_half[i].x <= list[end - 1].x))
					{
						by_y[count].x = other_half[i].x;
						by_y[count].y = other_half[i].y;
						count++;
					}
				}
			}
		}
		// print_list(by_y, length, false);

		// print_list(by_y, length, false);
		//

		// selecting root
		int m = median_index(length);
		location *rt = malloc(sizeof(location));
		if (gen % 2 == 0)
		{
			*rt = by_x[m];
		}
		else
		{
			*rt = by_y[m];
		}
		// printf("selecting (%lf, %lf) as the root\n\n", rt->x, rt->y);

		location *temp = loc_add(out, rt);
		out = temp;
		free(rt);

		if (length > 1)
		{
			if (gen % 2 == 0)
			{
				balance_locations(by_x, by_y, m + 1, length, length, length - m - 1, gen + 1, out, true);
				balance_locations(by_x, by_y, 0, m, length, m, gen + 1, out, true);

				free(by_x);
				free(by_y);
			}
			else
			{
				balance_locations(by_y, by_x, m + 1, length, length, length - m - 1, gen + 1, out, false);
				balance_locations(by_y, by_x, 0, m, length, m, gen + 1, out, false);

				free(by_y);
				free(by_x);
			}
		}
		else
		{
			free(by_x);
			free(by_y);
		}

		return out;
	}
	return NULL;
}

/**
 * Adds a copy of the given point to the given k-d tree.  There is no
 * effect if the point is already in the tree.  The tree need not be
 * balanced after the add.  The return value is true if the point was
 * added successfully and false otherwise (if the point was already in the
 * tree or if there was a memory allocation error when trying to add).
 *
 * @param t a pointer to a valid k-d tree, non-NULL
 * @param p a pointer to a valid location, non-NULL
 * @return true if and only if the point was successfully added
 */
bool kdtree_add(kdtree *t, const location *p)
{
	if (t->head == NULL)
	{
		location *p_copy = malloc(sizeof(location));
		p_copy->x = p->x;
		p_copy->y = p->y;

		node *rt = malloc(sizeof(node));
		// fprintf(stderr, "root node %p\n", rt);
		t->head = rt;

		t->head->root = p_copy;
		t->head->left = NULL;
		t->head->right = NULL;
		t->head->region = malloc(4 * sizeof(double));

		// free(p_copy);

		return true;
	}
	else
	{
		int add_aux = kdtree_add_aux(t->head, p, t->depth, 0);
		if (add_aux > t->depth)
		{
			t->depth = add_aux;
		}
		return (add_aux >= 0);
	}
}

// bool is_empty(location *loc)
// {
// 	return (loc->x == 0 && loc->y == 0);
// }

// returns depth of inserted node
int kdtree_add_aux(node *head, const location *p, int tree_depth, int count)
{
	location *p_copy = malloc(sizeof(location));
	p_copy->x = p->x;
	p_copy->y = p->y;

	if (loc_equals(head->root, p_copy))
	{
		free(p_copy);
		// printf("Location failed to add: %lf %lf\n", head->root->x, head->root->y);
		// printf("%d\n", -1);
		return -1;
	}

	node *new_node = malloc(sizeof(node));
	// fprintf(stderr, "allocing %p\n", new_node);
	new_node->left = NULL;
	new_node->right = NULL;
	new_node->region = malloc(4*sizeof(double));
	new_node->root = p_copy;

	if (count % 2 == 0)
	{
		if (p_copy->x < head->root->x)
		{	
			if (head->left != NULL)
			{
				free(new_node->region);
				free(new_node);
				free(p_copy);
				return kdtree_add_aux(head->left, p, tree_depth, count + 1);
			}
			else
			{
				head->left = new_node;
				// fprintf(stderr, "adding new_node to head\n");
			}
		}
		else
		{
			if (head->right != NULL)
			{
				free(new_node->region);
				free(new_node);
				free(p_copy);
				return kdtree_add_aux(head->right, p, tree_depth, count + 1);
			}
			else
			{
				head->right = new_node;
				// fprintf(stderr, "adding new_node to head\n");
			}
		}
	}
	else
	{
		if (p_copy->y < head->root->y)
		{
			if (head->left != NULL)
			{
				free(new_node->region);
				free(new_node);
				free(p_copy);
				return kdtree_add_aux(head->left, p, tree_depth, count + 1);
			}
			else
			{
				head->left = new_node;
				// fprintf(stderr, "adding new_node to head\n");
			}
		}
		else
		{
			if (head->right != NULL)
			{
				free(new_node->region);
				free(new_node);
				free(p_copy);			
				return kdtree_add_aux(head->right, p, tree_depth, count + 1);
			}
			else
			{
				head->right = new_node;
				// fprintf(stderr, "adding new_node to head\n");
			}
		}
	}

	return count + 1;
}

bool loc_equals(location *l1, location *l2)
{
	return (l1->x == l2->x && l1->y == l2->y);
}




/**
 * Determines if the given tree contains a point with the same coordinates
 * as the given point.
 *
 * @param t a pointer to a valid k-d tree, non-NULL
 * @param p a pointer to a valid location, non-NULL
 * @return true if and only of the tree contains the location
 */

bool kdtree_contains(const kdtree *t, const location *p)
{
	return kdtree_contains_aux(t->head, p);
}

bool kdtree_contains_aux(node *head, const location *p)
{
	location *p_copy = malloc(sizeof(location));
	p_copy->x = p->x;
	p_copy->y = p->y;

	if (head == NULL)
	{
		free(p_copy);
		return false;
	}
	else if (loc_equals(head->root, p_copy))
	{
		free(p_copy);
		return true;
	}
	else
	{
		free(p_copy);		
		return (kdtree_contains_aux(head->left, p) || kdtree_contains_aux(head->right, p));
	}
}

/**
 * Copies p's nearest neighbor to neighbor.  
 * Ties are broken arbitrarily.
 * There is no change to neighbor and distance is set to infinity if
 * the tree is empty.  If p is equal to a point in the tree then p is
 * copied into neighbor and distance is set to zero.
 *
 * @param t a pointer to a valid k-d tree, non-NULL
 * @param p a pointer to a valid location, non-NULL
 * @param d a pointer to a double, non-NULL
 * @param neighbor a pointer to a valid location, non-NULL
 */

void kdtree_nearest_neighbor(kdtree *t, const location *p, location *neighbor, double *d)
{
	// printf("Startup!\n");
	if (t->head == NULL)
	{
		*d = INFINITY;
		return;
	}

	location *p_copy = malloc(sizeof(location));
	p_copy->x = p->x;
	p_copy->y = p->y;

	produce_regions(t->head, 0, 0, 0, 0, 0, true);
	// printf("Regions produced\n");

	if (kdtree_contains(t, p))
	{
		neighbor->x = p_copy->x;
		neighbor->y = p_copy->y;
		*d = 0.0;
		free(p_copy);
		return;
	}

	// set starting neighbor as the leaf where p would have gone if not for it
	// printf("Neighbor!: (%lf, %lf)\n", neighbor->x, neighbor->y);
	location *start = starting_point(t->head, p_copy, 0, t->depth)->root;

	// printf("starting point determined: %lf %lf\n", start->x, start->y);


	// start->x = 10.0;
	// start->y = 1.0;
	double curr_min_dist = location_distance(p_copy, start);

	node *out_node = malloc(sizeof(node));
	out_node->root = NULL;
	out_node->left = NULL;
	out_node->right = NULL;
	out_node->region = NULL;

	out_node = find_neighbor(t->head, p_copy, start, &curr_min_dist, true, out_node);
	neighbor->x = out_node->root->x;
	neighbor->y = out_node->root->y;
	// printf("Neighbor: %lf %lf\n", neighbor->x, neighbor->y);
	*d = location_distance(p, neighbor);
	free(p_copy);
	free(out_node);
}

node *find_neighbor(node *head, location *p, location *curr_neighbor, double *curr_min_dist, bool is_head, node *out)
{
	// printf("Curr_min_i = %lf\n", *curr_min_dist);
	location *BL = malloc(sizeof(location));
	location *TR = malloc(sizeof(location));

	if (out->root == NULL)
	{
		out->root = curr_neighbor;
	}

	if (head != NULL)
	{
		BL->x = head->region[0];
		BL->y = head->region[1];
		TR->x = head->region[2];
		TR->y = head->region[3];
	}
	else
	{
		return NULL;
	}

	// printf("Regions of candidate (%lf %lf): %lf %lf %lf %lf\n", head->root->x, head->root->y, BL->x, BL->y, TR->x, TR->y);

	if (location_distance_to_rectangle(p, BL, TR) < *curr_min_dist || is_head)
	{
		// printf("Comparing point (%lf %lf) to potential neighbor (%lf %lf):\n", p->x, p->y, head->root->x, head->root->y);
		// printf("Distance: %lf\n", location_distance(head->root, p));
		if (location_distance(head->root, p) < *curr_min_dist)
		{
			*curr_min_dist = location_distance(head->root, p);
			curr_neighbor = head->root;
			out->root = curr_neighbor;
			// printf("Curr_min_f = %lf\n", *curr_min_dist);
		}

		if (head->right != NULL)
		{
			find_neighbor(head->right, p, curr_neighbor, curr_min_dist, false, out);
		}

		if (head->left != NULL)
		{
			find_neighbor(head->left, p, curr_neighbor, curr_min_dist, false, out);
		}
	}

	free(BL);
	free(TR);
	return out;
}

node *starting_point(node *head, location *p, int count, int depth)
{
	if (count < depth)
	{
		if (count % 2 == 0)
		{
			if (p->x < head->root->x)
			{
				if (head->left != NULL)
				{
					return starting_point(head->left, p, count + 1, depth);
				}
				else
				{
					// printf("Starting point: (%lf, %lf)\n", head->root->x, head->root->y);
					return head;
				}
			}
			else
			{
				if (head->right != NULL)
				{
					return starting_point(head->right, p, count + 1, depth);
				}
				else
				{
					// printf("Starting point: (%lf, %lf)\n", head->root->x, head->root->y);
					return head;
				}
			}
		}
		else
		{
			if (p->y < head->root->y)
			{
				if (head->left != NULL)
				{
					return starting_point(head->left, p, count + 1, depth);
				}
				else
				{
					// printf("Starting point: (%lf, %lf)\n", head->root->x, head->root->y);
					return head;
				}
			}
			else
			{
				if (head->right != NULL)
				{
					return starting_point(head->right, p, count + 1, depth);
				}
				else
				{
					// printf("Starting point: (%lf, %lf)\n", head->root->x, head->root->y);
					return head;
				}
			}
		}
	}
	else
	{
		// printf("Starting point: (%lf, %lf)\n", head->root->x, head->root->y);
		return head;
	}
}

void produce_regions(node *head, double BLx, double BLy, double TRx, double TRy, int gen, bool is_head)
{
	double A = head->root->x;
	double B = head->root->y;

	if (is_head)
	{
		BLx = -INFINITY;
		BLy = -INFINITY;
		TRx = INFINITY;
		TRy = INFINITY;

		head->region[0] = BLx;
		head->region[1] = BLy;
		head->region[2] = TRx;
		head->region[3] = TRy;
		// print_region(A, B, BLx, BLy, TRx, TRy, gen);
	}
	
	if (gen % 2 == 0)
	{
		if (head->right != NULL)
		{
			head->right->region[0] = A;
			head->right->region[1] = BLy;
			head->right->region[2] = TRx;
			head->right->region[3] = TRy;
			// location *right = head->right->root;
			// print_region(right->x, right->y, A, BLy, TRx, TRy, gen + 1);
			produce_regions(head->right, A, BLy, TRx, TRy, gen + 1, false);
		}
		if (head->left != NULL)
		{
			head->left->region[0] = BLx;
			head->left->region[1] = BLy;
			head->left->region[2] = A;
			head->left->region[3] = TRy;
			// location *left = head->left->root;
			// print_region(left->x, left->y, BLx, BLy, A, TRy, gen + 1);
			produce_regions(head->left, BLx, BLy, A, TRy, gen + 1, false);
		}
	}
	else
	{
		if (head->right != NULL)
		{
			head->right->region[0] = BLx;
			head->right->region[1] = B;
			head->right->region[2] = TRx;
			head->right->region[3] = TRy;
			// location *right = head->right->root;
			// print_region(right->x, right->y, BLx, B, TRx, TRy, gen + 1);
			produce_regions(head->right, BLx, B, TRx, TRy, gen + 1, false);
		}
		if (head->left != NULL)
		{
			head->left->region[0] = BLx;
			head->left->region[1] = BLy;
			head->left->region[2] = TRx;
			head->left->region[3] = B;
			// location *left =	 head->left->root;
			// print_region(left->x, left->y, BLx, BLy, TRx, B, gen + 1);
			produce_regions(head->left, BLx, BLy, TRx, B, gen + 1, false);
		}
	}
}

void print_region(double A, double B, double BLx, double BLy, double TRx, double TRy, int gen)
{
	printf("%d (%.3lf, %.3lf) [(%.3lf, %.3lf)-(%.3lf, %.3lf)]\n", gen, A, B, BLx, BLy, TRx, TRy);
}

/**
 * Passes the points in the given tree to the given function
 * in an arbitrary order.  The last argument to this function is also passed
 * to the given function along with each point.
 *
 * @param t a pointer to a valid k-d tree, non-NULL
 * @param f a pointer to a function that takes a location and an extra
 * argument, and does not modify t, non-NULL
 * @param arg a pointer
 */
void kdtree_for_each(const kdtree* r, void (*f)(const location *, void *), void *arg)
{
	traverse(r->head, f, arg);
}

void traverse(node *head, void (*f)(const location *, void *), void *arg)
{
	if (head != NULL)
	{
		traverse(head->right, f, arg);
		f(head->root, arg);
		traverse(head->left, f, arg);
	}
}


void print_root(const location *root, void *arg)
{
	printf("(%lf, %lf)\n", root->x, root->y);
}

void print_kdtree(kdtree *t)
{
	kdtree_for_each(t, print_root, NULL);
}

/**
 * Destroys the given k-d tree.  The tree is invalid after being destroyed.
 *
 * @param t a pointer to a valid k-d tree, non-NULL
 */
void kdtree_destroy(kdtree *t)
{
	traverse_free(t->head);
	free(t);
}

void traverse_free(node *head)
{
	if (head == NULL)
	{
		return;
	}
	if (head->left != NULL)
	{
		traverse_free(head->left);
	}
	if (head->right != NULL)
	{
		traverse_free(head->right);
	}
	free(head->root);
	free(head->region);
	// fprintf(stderr, "freeing %p\n", head);
	free(head);	
}