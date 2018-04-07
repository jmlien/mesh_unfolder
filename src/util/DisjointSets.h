#pragma once

#include <vector>
#include <iostream>
#include <cassert>
using namespace std;
// DisjointSets class
//
// CONSTRUCTION: with int representing initial number of sets
//
// ******************PUBLIC OPERATIONS*********************
// void union( root1, root2 ) --> Merge two sets
// int find( x )              --> Return set containing x
// ******************ERRORS********************************
// Error checking or parameters is performed


/**
 * Disjoint set class, using union by rank
 * and path compression.
 * Elements in the set are numbered starting at 0.
 * @author Mark Allen Weiss (modified to C++ by JML)
 */

class DisjointSets
{
public:

     DisjointSets( int size )
     {
        this->n_of_sets = this->cap = size;
        this->s = new int [ size ];
        assert(this->s);
        for( int i = 0; i < n_of_sets; i++ ) s[ i ] = -1;
    }

    /**
     * Union two disjoint sets using the height heuristic.
     * root1 and root2 are distinct and represent set names.
     * @param root1 the root of set 1.
     * @param root2 the root of set 2.
     * @throws IllegalArgumentException if root1 or root2
     * are not distinct roots.
     */
    int unite( int root1, int root2 )
    {
        assertIsRoot( root1 );
        assertIsRoot( root2 );
        if( root1 == root2 )
        {
          cerr<< "! Error: DisjointSets: Union: root1 == root2: " << root1<<endl;
          assert(false);
          return -1;
        }

        n_of_sets--; //decrease the number of sets

        if( s[ root2 ] <= s[ root1 ] )  // root2 is larger
        {
            if(s[ root2 ] == s[ root1 ]) s[ root2 ]--; // Update size
            s[ root1 ] = root2;     // Make root2 new root
            return root2;
        }
        else
        {
            s[ root2 ] = root1;        // Make root1 new root
            return root1;
        }
    }

    /**
     * Perform a find with path compression.
     * @param x the element being searched for.
     * @return the set containing x.
     * @throws IllegalArgumentException if x is not valid.
     */
    int find( int x )
    {
        assertIsItem( x );
        if( s[ x ] < 0 )
            return x;
        else
            return s[ x ] = find( s[ x ] );
    }

    // return the number of disjoint sets remaining
    int getNumSets()
    {
      return n_of_sets;
    }

private:

    int * s;
    int n_of_sets, cap;

    void assertIsRoot( int root )
    {
        assertIsItem( root );
        if( s[ root ] >= 0 )
        {
          cerr<<"! Error: DisjointSets: "<<  root << " is not a root, its parent is: "<<s[ root ]<<endl;
          assert(false);
        }
    }

    void assertIsItem( int x )
    {
        if( x < 0 || x >= cap )
        {
            cerr<<"! Error: DisjointSets: "<<  x<<" is not an item"<<endl;
        }
    }

};
