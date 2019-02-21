#include <iostream>    
#include <stdio.h>

#include "hycon3d.h"

#ifndef _MYBOXLIST_H
#define _MYBOXLIST_H 1

/*
int cmppointx(void* a,void* b){
	return (int)(((point_avl*)a)->x - ((point_avl*)b)->x);
}
*/

class double_linked
{
 public:
    struct node
    {
        box data;
        node* prev;
        node* next;
        node(box t, node* p, node* n) : data(t), prev(p), next(n) {}
    };
    node* head;
    node* tail;

    double_linked() : head( NULL ), tail ( NULL ) {}


    bool empty() const { return ( !head || !tail ); }
    operator bool() const { return !empty(); } 
    void push_back(box);
    void push_front(box);
    box pop_back();
    box pop_front();
    void print_boxes();

    ~double_linked()
    {
        while(head)
        {
            node* temp(head);
            head=head->next;
            delete temp;
        }
    }
};

void printbox(box b);

#endif
