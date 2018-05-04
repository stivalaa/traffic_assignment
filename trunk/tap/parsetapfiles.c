/*****************************************************************************
 * 
 * File:    parsetapfiles.c
 * Author:  Alex Stivala
 * Created: March 2011
 *
 * $Id: parsetapfiles.c 696 2011-09-14 02:02:22Z astivala $
 *
 * Functions to parse the test data _trips, _net, (traffic assignment input)
 * data in the format from
 *
 * http://www.bgu.ac.il/~bargera/tntp/
 *
 ****************************************************************************/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "parsetapfiles.h"


#define MAX(a,b) ((a) > (b) ? (a) : (b))

/****************************************************************************
 *
 * constants and type definitions
 *
 ****************************************************************************/

#define MAX_LINE_LEN 4096 /* size of line buffer */

/****************************************************************************
 *
 * externally visible functions
 *
 ****************************************************************************/

/* 
 *
 * parse_net_file() - 
 *   Parse the network file as described at http://www.bgu.ac.il/~bargera/tntp/
 *
 *   E.g.
 *   <NUMBER OF ZONES> 24
 *   <NUMBER OF NODES> 24
 *   <FIRST THRU NODE> 1
 *   <NUMBER OF LINKS> 76
 *   <END OF METADATA>
 *   
 *   ~ 	Init node 	Term node 	Capacity 	Length 	Free Flow Time 	B *	Power	Speed limit 	Toll 	Type	;
 *	1	2	25900.20064	6	6	0.15	4	0	0
 *
 *   Parameters:
 *      fp - open (read) filename to parse from
 *      net (OUT) - structure describing network parsed form file, inc Link list
 *
 *   Return value:
 *       0 if Ok nonzero for error
 *
 */
int parse_net_file(FILE *fp, net_data_t *net)
{
  long link_count = 0;
  char buf[MAX_LINE_LEN];
  long max_node_num = 0;
  const int INITIAL_ALLOC = 1000; /* alloc this any links initially */
  static const char *zones_tag = "<NUMBER OF ZONES>";
  static const char *nodes_tag = "<NUMBER OF NODES>";
  static const char *thru_tag  = "<FIRST THRU NODE>";
  static const char *links_tag = "<NUMBER OF LINKS>";
  static const char *end_tag  =  "<END OF METADATA>";

  if (!(net->links = (link_data_t *)malloc(INITIAL_ALLOC * sizeof(link_data_t))))
  {
    fprintf(stderr, "malloc net->links failed\n");
    return -1;
  }

  while (fgets(buf, sizeof(buf)-1, fp))
  {
    if (buf[0]  == '~' || strspn(buf, " \t\n\r") == strlen(buf))
      continue; /* skip comment or whitespace lines */
    if (buf[0] == '<') /* metadata */
    {
      char *closeangle = strchr(buf, '>');
      char *val;
      if (!closeangle)
      {
        fprintf(stderr, "error: no closing '>' on metadata tag: %s\n",buf);
        return -1;
      }
      val = closeangle+1 + strcspn(closeangle+1, " \t");
      if (strncmp(buf, zones_tag, strlen(zones_tag)) == 0)
        net->num_zones = atol(val);
      else if (strncmp(buf, nodes_tag, strlen(nodes_tag)) == 0)
        net->num_nodes = atol(val);
      else if (strncmp(buf, thru_tag, strlen(thru_tag)) == 0)
        net->first_thru_node = atol(val);
      else if (strncmp(buf, links_tag, strlen(links_tag)) == 0)
        net->num_links = atol(val);
      else if (strncmp(buf, end_tag, strlen(end_tag)) == 0)
        continue;
      else
      {
        fprintf(stderr, "warning: unknown metadata tag: %s\n", buf);
      }
    }
    else
    {
      link_data_t *link;
      char *saveptr = NULL;
      if (link_count >= INITIAL_ALLOC)
      {
        if (!(net->links = (link_data_t *)realloc(net->links, 
                                          (link_count+1)*sizeof(link_data_t))))
        {
          fprintf(stderr, "realloc net->links failed\n");
          return -1;
        }
      }
      link = &net->links[link_count];
      link->init_node = atol(strtok_r(buf, " \t", &saveptr));
      max_node_num = MAX(max_node_num, link->init_node);
      link->term_node = atol(strtok_r(NULL, " \t", &saveptr));
      max_node_num = MAX(max_node_num, link->term_node);
      link->capacity = atof(strtok_r(NULL, " \t", &saveptr));
      link->length = atof(strtok_r(NULL, " \t", &saveptr));
      link->free_flow_time = atof(strtok_r(NULL, " \t", &saveptr));
      link->B = atof(strtok_r(NULL, " \t", &saveptr));
      link->power = atof(strtok_r(NULL, " \t", &saveptr));
      link->speed_limit = atof(strtok_r(NULL, " \t", &saveptr));
      link->toll = atof(strtok_r(NULL, " \t", &saveptr));
      link->link_type = (short)atoi(strtok_r(NULL, " \t", &saveptr));
      link_count++;
    }
  }

  if (link_count != net->num_links)
  {
    fprintf(stderr, "error: num_links %ld but %ld links\n",
            net->num_links, link_count);
    return -1;
  }
  return 0;
}



/*
 * parse_trips_file() - parse O-D demand (trips) file
 *
 *    Parse the origin-destination demand file (trips file) which gives
 *   the demand fro travel between origin and destination (zone) pairs
 *   as described at http://www.bgu.ac.il/~bargera/tnt/p
 *
 *   Paraemeters:
 *      fp - open (read) filepointer to read from
 *      demands - (OUT) array of array of dmeand data structs parsed from file
 *                for each origin demands[origin] is an array of demands structs
 *                terminated by one with 0  dest (origins and dests start at 1)
 *      num_zones - (OUT) number of origin entries in demands array
 *
 *   Return value:
 *      0 if OK else nonzero
 *
 *
 *      
 *   Format of file  is e.g.:
 *
 *   <NUMBER OF ZONES> 24
 *   <TOTAL OD FLOW> 360600.0
 *   <END OF METADATA>
 *
 *
 *   Origin  1
 *       1 :      0.0;     2 :    100.0;     3 :    100.0;     4 :    500.0;     5 :    200.0;
 *       6 :    300.0;     7 :    500.0;     8 :    800.0;     9 :    500.0;    10 :   1300.0;
 *
 */
int parse_trips_file(FILE *fp, demand_data_t **demands[], long *num_zones)

{
  const int INITIAL_ALLOC = 10;
  long demand_count = 0;
  long origin_demand_count = 0;
  double total_demands = 0;
  char buf[MAX_LINE_LEN];
  double total_od_flow = 0;
  long cur_origin = 0;
  static const char *zones_tag = "<NUMBER OF ZONES>";
  static const char *total_tag = "<TOTAL OD FLOW>";
  static const char *end_tag  =  "<END OF METADATA>";

  *demands = NULL;
  while (fgets(buf, sizeof(buf)-1, fp))
  {
    if (buf[0]  == '~' || strspn(buf, " \t\n\r") == strlen(buf))
      continue; /* skip comment or whitespace lines */
    if (buf[0] == '<') /* metadata */
    {
      char *closeangle = strchr(buf, '>');
      char *val;
      if (!closeangle)
      {
        fprintf(stderr, "error: no closing '>' on metadata tag: %s\n",buf);
        return -1;
      }
      val = closeangle+1 + strcspn(closeangle+1, " \t");
      if (strncmp(buf, zones_tag, strlen(zones_tag)) == 0)
        *num_zones = atol(val);
      else if (strncmp(buf, total_tag, strlen(total_tag)) == 0)
        total_od_flow = atof(val);
      else if (strncmp(buf, end_tag, strlen(end_tag)) == 0)
      {
        if (!*demands)
        {
          if (!(*demands = (demand_data_t **)calloc(*num_zones + 1,
                                                    sizeof(demand_data_t*))))
          {
            fprintf(stderr, "calloc demands failed\n");
            return -1;
          }
        }
      }
      else
      {
        fprintf(stderr, "warning: unknown metadata tag: %s\n", buf);
      }
    }
    else if (strncmp(buf, "Origin", strlen("Origin")) == 0)
    {
      if (sscanf(buf, "Origin %ld\n", &cur_origin) != 1)
      {
        fprintf(stderr, "error parsing origin: %s\n", buf);
        return -1;
      }
      if (cur_origin > *num_zones)
        fprintf(stderr, "warning: origin %ld > num_zones %ld\n",
                cur_origin, *num_zones);
      if ((*demands)[cur_origin])
      {
        fprintf(stderr, "error: multiple entries for Origin %ld\n", cur_origin);
        return -1;
      }
      if (!((*demands)[cur_origin] = calloc(INITIAL_ALLOC, sizeof(demand_data_t))))
      {
        fprintf(stderr, "calloc for origin %ld failed\n", cur_origin);
        return -1;
      }
      origin_demand_count = 0; /* entries for cur_origin */
    }
    else
    {
      demand_data_t *demand;
      char *saveptr1 = NULL;
      char *saveptr2 = NULL;
      char *str1, *token;
      char *dest_str;
      for (str1 = buf; /* nothing */; str1 = NULL) 
      {
        if (origin_demand_count >= INITIAL_ALLOC)
        {
          if (!((*demands)[cur_origin] = (demand_data_t *)
                realloc((*demands)[cur_origin],
                        (origin_demand_count+1)*sizeof(demand_data_t))))
          {
            fprintf(stderr, "realloc demands failed\n");
            return -1;
          }
          ((*demands)[cur_origin][origin_demand_count]).dest = 0;/* end marker */
        }
        demand = &((*demands)[cur_origin][origin_demand_count]);
        token = strtok_r(str1, ";", &saveptr1);
        if (token == NULL)
          break;
        dest_str = strtok_r(token, ":", &saveptr2);
        if (!dest_str || strchr(dest_str, '\n') != NULL)
          continue;
        demand->dest = atol(dest_str);
        if (demand->dest > *num_zones)
          fprintf(stderr, "warning: dest %ld > num_zones %ld\n",
                  demand->dest, *num_zones);
        demand->demand = atof(strtok_r(NULL,  ":", &saveptr2));
        origin_demand_count++;
        total_demands += demand->demand;
        demand_count++;
      }
    }
  }

  /* sometimnes there are missing origins e.g. in Philadelphia data */
  /* we need to allocate just the 0 terminator so they don't cause a problem*/
  for (cur_origin = 1; cur_origin <= *num_zones; cur_origin++)
  {
    if (!(*demands)[cur_origin])
    {
      fprintf(stderr, "warning: no trip data for origin %ld\n", cur_origin);
      if (!((*demands)[cur_origin] = calloc(1, sizeof(demand_data_t))))
      {
        fprintf(stderr, "calloc for origin %ld failed\n", cur_origin);
        return -1;
      }
    }
  }
 
  if (fabs(total_od_flow - total_demands) > 1e-08)
    fprintf(stderr, "warning: total_od_flow %f != %f\n", 
            total_od_flow, total_demands);
  return 0;
}


/*
 * parse_flows_file() - Parse the traffic assignment output flows file
 *
 *  """
 *   Parse the flow file giving volume and cost from output from traffic assign
 *   
 *   Paraemters:
 *      flow_fp - open (read) filehandle of the
 *                 flow file output by traffic assignment program
 *      link_flows (OUT) array of link_flow_t structures allocated and
 *                 parsed from flow_fp file
 *      num_links (OUT) number of elements in link_flows array
 *      
 *
 *   Return value:
 *      0 if OK else nonzero on error
 *
 * File has format e.g.:
 *
 *      <NUMBER OF NODES>       25
 *      <NUMBER OF LINKS>       76
 *      <END OF METADATA>
 *
 *
 *      ~       Tail    Head    :       Volume  Cost    ;
 *              1       2       :       4494.5008499645437041   6.0008161234622576785   ;
 *               1       3       :       8119.1900669362912595   4.0086912217072878661  ;
 *
 */

int parse_flows_file(FILE *flow_fp, link_flow_t **link_flows, long *num_links)
{
  long link_count = 0;
  char buf[MAX_LINE_LEN];
  long max_node_num = 0;
  long meta_num_nodes = 0, meta_num_links =0;
  const int INITIAL_ALLOC = 1000; /* alloc this any links initially */
  static const char *nodes_tag = "<NUMBER OF NODES>";
  static const char *links_tag = "<NUMBER OF LINKS>";
  static const char *end_tag  =  "<END OF METADATA>";

  if (!(*link_flows = (link_flow_t *)malloc(INITIAL_ALLOC*sizeof(link_flow_t))))
  {
    fprintf(stderr, "malloc link_flows failed\n");
    return -1;
  }

  while (fgets(buf, sizeof(buf)-1, flow_fp))
  {
    if (buf[0]  == '~' || strspn(buf, " \t\n\r") == strlen(buf))
      continue; /* skip comment or whitespace lines */
    if (buf[0] == '<') /* metadata */
    {
      char *closeangle = strchr(buf, '>');
      char *val;
      if (!closeangle)
      {
        fprintf(stderr, "error: no closing '>' on metadata tag: %s\n",buf);
        return -1;
      }
      val = closeangle+1 + strcspn(closeangle+1, " \t");
      if (strncmp(buf, nodes_tag, strlen(nodes_tag)) == 0)
        meta_num_nodes = atol(val);
      else if (strncmp(buf, links_tag, strlen(links_tag)) == 0)
        meta_num_links = atol(val);
      else if (strncmp(buf, end_tag, strlen(end_tag)) == 0)
        continue;
      else
      {
        fprintf(stderr, "warning: unknown metadata tag: %s\n", buf);
      }
    }
    else
    {
      link_flow_t *link;
      char *saveptr = NULL;
      if (link_count >= INITIAL_ALLOC)
      {
        if (!(*link_flows = (link_flow_t *)realloc(*link_flows, 
                                          (link_count+1)*sizeof(link_flow_t))))
        {
          fprintf(stderr, "realloc link_flows failed\n");
          return -1;
        }
      }
      link = &((*link_flows)[link_count]);
      link->init_node = atol(strtok_r(buf, "\t", &saveptr));
      max_node_num = MAX(max_node_num, link->init_node);
      link->term_node = atol(strtok_r(NULL, "\t", &saveptr));
      max_node_num = MAX(max_node_num, link->term_node);
      if (strchr(strtok_r(NULL, "\t", &saveptr), ":") != 0)
      {
        fprintf(stderr, "bad line in flows file: %s\n", buf);
        return -1;
      }
      link->volume = atof(strtok_r(NULL, "\t", &saveptr));
      link->cost = atof(strtok_r(NULL, "\t", &saveptr));
      link_count++;
    }
  }

  if (link_count != meta_num_links)
  {
    fprintf(stderr, "error: num_links %ld but %ld links\n",
            meta_num_links, link_count);
    return -1;
  }
  *num_links = link_count;
  return 0;
}

/* 
 *
 * parse_net_mod_file() - 
 *   Parse the network modifications file, which describes changes to
 *   the net file.
 *
 *   This has the format (one change per line):
 *   
 *    ChangeId ChangeType Initnode  Termnode  Capacity 	Length 	FreeFlowTime  B Power	Speedlimit   Toll Type	
 *
 *   The ChangeId must be suitable to be put in a filename as flows output
 *    files have it appended. 
 *   Note that multiple entries with a single ChangeId are allowed:
 *   this descrbies multiple modifications considered as part of a single
 *   indivisible network improvement. This is often used for example
 *   to add bidirectional links for a new road.
 *
 *   The ChangeType can be 'ADD' or 'CHANGE'
 *
 *   Note (at least for now) nodes have to already exist (ie cannot add nodes,
 *   only links between existing nodes). TODO allow new nodes. 
 *   TODO allow link deletion.
 *
 *   Parameters:
 *      fp - open (read) filename to parse from
 *      mods (OUT) - array of structures describing network modifications
 *      num_mods (OUT) - number of entries in mods array
 *
 *   Return value:
 *      0 if OK else nonzero on error
 *
 */
int parse_net_mod_file(FILE *fp, net_mod_t **mods, int *num_mods)
{
  long mod_count = 0;
  char buf[MAX_LINE_LEN];
  const int INITIAL_ALLOC = 10; /* alloc this any mods initially */
  link_data_t *link;
  char *changetype_str;


  if (!(*mods = (net_mod_t *)malloc(INITIAL_ALLOC * sizeof(net_mod_t))))
  {
    fprintf(stderr, "malloc mods failed\n");
    return -1;
  }

  while (fgets(buf, sizeof(buf)-1, fp))
  {
    if (buf[0]  == '~' || strspn(buf, " \t\n\r") == strlen(buf))
      continue; /* skip comment or whitespace lines */

    net_mod_t *mod;
    char *saveptr = NULL;
    if (mod_count >= INITIAL_ALLOC)
    {
      if (!(*mods = (net_mod_t *)realloc(*mods,
                                         (mod_count+1)*sizeof(net_mod_t))))
      {
        fprintf(stderr, "realloc mods failed\n");
        return -1;
      }
    }
    mod = &((*mods)[mod_count]);
    strncpy(mod->change_id, strtok_r(buf, "\t", &saveptr), MAX_CHANGE_ID_LEN);
    changetype_str = strtok_r(NULL, "\t", &saveptr);
    if (strcmp(changetype_str, "ADD") == 0)
      mod->modtype = MOD_TYPE_ADD_E;
    else if (strcmp(changetype_str, "CHANGE") == 0)
      mod->modtype = MOD_TYPE_CHANGE_E;
    else
    {
      fprintf(stderr, "invalid change type '%s'\n", changetype_str);
      return -1;
    }
    link = &mod->mod_link;
    link->init_node = atol(strtok_r(NULL, "\t", &saveptr));
    link->term_node = atol(strtok_r(NULL, "\t", &saveptr));
    link->capacity = atof(strtok_r(NULL, "\t", &saveptr));
    link->length = atof(strtok_r(NULL, "\t", &saveptr));
    link->free_flow_time = atof(strtok_r(NULL, "\t", &saveptr));
    link->B = atof(strtok_r(NULL, "\t", &saveptr));
    link->power = atof(strtok_r(NULL, "\t", &saveptr));
    link->speed_limit = atof(strtok_r(NULL, "\t", &saveptr));
    link->toll = atof(strtok_r(NULL, "\t", &saveptr));
    link->link_type = (short)atoi(strtok_r(NULL, "\t", &saveptr));
    mod->project_cost = atof(strtok_r(NULL, "\t", &saveptr));
    ++mod_count;
  }
  *num_mods = mod_count;
  return 0;
}


/* 
 * link_data_compar() - qsort comparison function for link_data_t
 * 
 * Compares by 'from' node number first  then by 'to' node number if equal
 * NB this must match the linkflow_entry_compar() and adjlist_entry_compar()
 * - see the comment block at top of tap_frankwolfe_pthread.c
 * re keeping arrays of link data sorted in this way.
 *
 */
int link_data_compar(const void *ent1, const void *ent2)
{
  const link_data_t *e1 = (const link_data_t *)ent1;
  const link_data_t *e2 = (const link_data_t *)ent2;
  
  if (e1->init_node < e2->init_node)
    return -1;
  else if(e1->init_node > e2->init_node)
    return 1;
  else
    return ( e1->term_node < e2->term_node ? -1 : 
             (e1->term_node > e2->term_node ? 1 : 0) );
}
  

/* 
 *
 * parse_trip_mod_file() - 
 *   Parse the O-D demand data modifications file, which describes changes to
 *   the trips file.
 *
 *   This has the format (one change per line):
 *   
 *    ChangeId ChangeType Origin Dest Multiplier
 *
 *   Note that multiple entries with a single ChangeId are allowed.
 *   The ChangeType can be 'P2P' for point-to-point (origin-dest) change
 *   or 'UNIFORM' for one-to-many or many-to-one change,in which case
 *   the multiplier is applied to all trips from orig (or to dest).
 *   There can also be 'ALL' in which case both orig and dest are 0 and
 *   all O-D demand data is multiplid by the multiplier value.
 *   Origin is th etrip origin node (zone) number, or 0 for uniform dest.
 *   Dest is the trip destination node (zone) number, or 0 for uniform orig.
 *   Multiplier is the deamdn multiplier e.g. 1.50 for 50% increase.
 *
 *
 *   Parameters:
 *      fp - open (read) filename to parse from
 *      mods (OUT) - array of structures describing network modifications
 *      num_mods (OUT) - number of entries in mods array
 *
 *   Return value:
 *      0 if OK else nonzero on error
 *
 */
int parse_trip_mod_file(FILE *fp, trip_mod_t **mods, int *num_mods)
{
  long mod_count = 0;
  char buf[MAX_LINE_LEN];
  const int INITIAL_ALLOC = 10; /* alloc this any mods initially */
  char *changetype_str;


  if (!(*mods = (trip_mod_t *)malloc(INITIAL_ALLOC * sizeof(trip_mod_t))))
  {
    fprintf(stderr, "malloc trip mods failed\n");
    return -1;
  }

  while (fgets(buf, sizeof(buf)-1, fp))
  {
    if (buf[0]  == '~' || strspn(buf, " \t\n\r") == strlen(buf))
      continue; /* skip comment or whitespace lines */

    trip_mod_t *mod;
    char *saveptr = NULL;
    if (mod_count >= INITIAL_ALLOC)
    {
      if (!(*mods = (trip_mod_t *)realloc(*mods,
                                         (mod_count+1)*sizeof(trip_mod_t))))
      {
        fprintf(stderr, "realloc trip mods failed\n");
        return -1;
      }
    }
    mod = &((*mods)[mod_count]);
    strncpy(mod->trip_change_id, strtok_r(buf,"\t",&saveptr), MAX_CHANGE_ID_LEN);
    changetype_str = strtok_r(NULL, "\t", &saveptr);
    if (strcmp(changetype_str, "P2P") == 0)
      mod->trip_modtype = TRIP_MOD_TYPE_POINT_E;
    else if (strcmp(changetype_str, "UNIFORM") == 0)
      mod->trip_modtype = TRIP_MOD_TYPE_UNIFORM_E;
    else if (strcmp(changetype_str, "ALL") == 0)
      mod->trip_modtype = TRIP_MOD_TYPE_ALL_E;
    else
    {
      fprintf(stderr, "invalid trip data change type '%s'\n", changetype_str);
      return -1;
    }
    mod->origin = atol(strtok_r(NULL, "\t", &saveptr));
    mod->dest = atol(strtok_r(NULL, "\t", &saveptr));
    mod->multiplier = atof(strtok_r(NULL, "\t", &saveptr));
    if (mod->trip_modtype == TRIP_MOD_TYPE_POINT_E && 
        (mod->origin == 0 || mod->dest == 0))
    {
      fprintf(stderr, 
              "ERROR: origin and dest must be nonzero for point-to-point\n");
      return -1;
    }
    if (mod->trip_modtype != TRIP_MOD_TYPE_ALL_E &&
        mod->origin == 0 && mod->dest == 0)
    {
      fprintf(stderr, "ERROR: cannot have 0 for both origin and dest except type ALL\n");
      return -1;
    }
    ++mod_count;
  }
  *num_mods = mod_count;
  return 0;

}
