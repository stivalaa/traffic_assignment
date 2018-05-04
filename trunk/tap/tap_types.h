#ifndef TAP_TYPES_H
#define TAP_TYPES_H
/*****************************************************************************
 * 
 * File:    tap_types.h
 * Author:  Alex Stivala
 * Created: March 2011
 *
 * $Id: tap_types.h 682 2011-09-12 06:40:52Z astivala $
 *
 * type defintions for TAP data (traffic assignment input)
 *
 ****************************************************************************/




/* represents each link from the net file. NB node numbers are numbered from 1
   BPR function: travel_time(Q) = T_0 * (1 + alpha*(Q/Q_max)^beta) */
typedef struct link_data_s
{
  long   init_node;          /* from node number */
  long   term_node;          /* to node number */
  double capacity;           /* Q_max in BPR function */
  double length;
  double free_flow_time;     /* T_0 in BPR function */
  double B;                  /* alpha in BPR function */
  double power;              /* beta in BPR function */
  double speed_limit;
  double toll;
  short  link_type;
} link_data_t;

/* represnets the network data frmo the net file */
typedef struct net_data_s
{
  /* metadata */
  long num_zones;            /*  number of zones */
  long num_nodes;            /* number of nodes */
  long first_thru_node;      /* first node that is not a zone (NB starts at 1) */
  long num_links;            /* number of links */
  /* link data */
  link_data_t *links;       /* array of num_links link structs */
} net_data_t;

/* represents a trip in the origin-destination demand (trips) file */
/* there will be an array of these for each origin to make iteration */
/* over origins and within them destionatinos efficient */
typedef struct demand_data_s
{
  long   dest;              /* destinatino node number (zone) (starts at 1) */
  double demand;            /* demand from orig to dest */
} demand_data_t;


/* represents flow data parsed from the traffic assingment output (flows) file */
typedef struct link_flow_s
{
  long   init_node;      /*  from node number */
  long   term_node;      /* to node number */
  double volume;         /* volume on link from TAP output */
  double cost;           /* cost on link from TAP output */
} link_flow_t;


#define MAX_CHANGE_ID_LEN 127 /* max length of a change identifier */

typedef enum net_mod_type_e
{
  MOD_TYPE_INVALID_E = 0,   /* invalid */
  MOD_TYPE_ADD_E     = 1,   /* add a link */
  MOD_TYPE_CHANGE_E  = 2    /* change link attribute(s) */
} net_mod_type_e;

/* represents a modification to the network, parsed from netmods file */
typedef struct net_mod_s
{
  char change_id[MAX_CHANGE_ID_LEN+1]; /* change identifier */
  net_mod_type_e modtype;  /* type of modification (add, change,...) */
  link_data_t mod_link;      /* change dor added link data */
  double project_cost;     /* cost of building this modification */
} net_mod_t;

typedef enum trip_mod_type_e
{
  TRIP_MOD_TYPE_INVALID_E = 0,  /* invalid */
  TRIP_MOD_TYPE_POINT_E   = 1,  /* point-to-point */
  TRIP_MOD_TYPE_UNIFORM_E = 2,  /* uniform: all dests for this origin change */
  TRIP_MOD_TYPE_ALL_E     = 3   /* all: all trip data multiplied by this */
} trip_mod_type_e;

/* represenets a chnage in origin-destination demand, parsed from tripmods file*/
typedef struct trip_mod_s
{
  char trip_change_id[MAX_CHANGE_ID_LEN+1]; /* trip modificatin identifier */
  trip_mod_type_e trip_modtype; /* type (uniform, point-to-point,...) */
  long origin;    /* trip origin node (zone) number or 0 for uniform dest*/
  long dest;      /* trip destination node (zone) number or 0 for uniform orig */
  double multiplier; /* multiply demand by this, e.g. 1.50 for 50% increase */
} trip_mod_t;

#endif /* TAP_TYPES_H */

