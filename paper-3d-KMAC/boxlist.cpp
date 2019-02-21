#include "boxlist.h"

void double_linked::push_back(box data)
{
    tail = new node(data, tail, NULL);
    if( tail->prev )
        tail->prev->next = tail;

    if( empty() )
        head = tail;
}

void double_linked::push_front(box data)
{
    head = new node(data, NULL, head);
    if( head->next )
        head->next->prev = head;

    if( empty() )
        tail = head;
}

box double_linked::pop_back()
{
    if( empty() )
        throw("double_linked : list empty");
    node* temp(tail);
    box data( tail->data );
    tail = tail->prev ;

    if( tail )
        tail->next = NULL;
    else
        head = NULL ;

    delete temp;
    return data;
}

box double_linked::pop_front()
{
    if( empty() )
        throw("double_linked : list empty");
    node* temp(head);
    box data( head->data );
    head = head->next ;

    if( head )
        head->prev = NULL;
    else
        tail = NULL;

    delete temp;
    return data;
}

void double_linked::print_boxes()
{
    if( empty() )
    {
        printf("double_linked : list empty");
    }
    else
        {
          node* temp(head);
          printf( "\n L=[...\n");
          while (temp)
          {
                   box data( temp->data );
                   printf("%f, %f, %f, %f, %f, %f", 
                      data.lx, data.ly, data.lz, 
                      data.ux, data.uy, data.uz);
                   temp = temp->next ;
                   if (temp) printf(";...\n");
          }
          delete temp;
          printf("]\n");
        }    
}

void printbox(box data)
{
      printf("lx=%f, ly=%f, lz=%f,\n ux=%f, uy=%f, uz=%f\n", 
                      data.lx, data.ly, data.lz, 
                      data.ux, data.uy, data.uz);
}
