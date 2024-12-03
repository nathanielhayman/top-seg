#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#include <cblas.h>

#define DIM 3

#ifdef DEBUG
/* When debugging is enabled, these form aliases to useful functions */
#define dbg_printf(...) printf(__VA_ARGS__)
#define dbg_requires(...) assert(__VA_ARGS__)
#define dbg_assert(...) assert(__VA_ARGS__)
#define dbg_ensures(...) assert(__VA_ARGS__)
#else
/* When debugging is disnabled, no code gets generated for these */
#define dbg_printf(...)
#define dbg_requires(...)
#define dbg_assert(...)
#define dbg_ensures(...)
#endif

struct path {
    int *points;
    double plen;
    int dlen;
};

typedef struct path path_t;

struct point {
    double *data;
    int seg_i;
};

typedef struct point pt_t;

struct segmented {
    pt_t *points;
    size_t rsize;
    int n_seg;
    path_t *paths;
};

typedef struct segmented seg_t;

struct mesh {
    const double **verts; 
    size_t vsize;
    const int* faces;
    size_t fsize;
};

typedef struct mesh mesh_t;

void print_seg(const seg_t *seg) {
    printf("============================\n");
    printf("Seglist at %d\n", seg);
    printf("Length: %d\n", seg->n_seg);
    printf("Num points: %d\n", seg->rsize);
    printf("Paths at %d\n", seg->paths);
    printf("First path length: %d\n", seg->paths != NULL ? seg->paths[0].dlen : 0);
    printf("============================\n");
}

void print_arr(const int *arr, int asize) {
    for (int i = 0; i < asize; i++) {
        printf("%d ", arr[i]);
    }
    printf("\n");
}

double vec_dist(const double *v1, const double *v2) {
    // Eigen::Vector3d x1(v1[0], v1[1], v1[2]);
    // Eigen::Vector3d x2(v2[0], v2[1], v2[2]);

    // return (x1 - x2).norm();

    double *im = (double *)malloc(DIM * sizeof(double));

    cblas_dcopy(DIM, v1, 1, im, 1);

    cblas_daxpy(DIM, -1.0, v2, 1, im, 1);

    double res = cblas_dnrm2(DIM, im, 1);

    free(im);

    return res;
}

void pt_cpy(pt_t *pt_dest, pt_t *pt_src) {
    int len = DIM * sizeof(double);
    pt_dest->data = malloc(len);
    memcpy(pt_dest->data, pt_src->data, len);
    pt_dest->seg_i = pt_src->seg_i;
}

void path_cpy(path_t *path_dest, path_t *path_src) {
    int len = path_src->dlen * sizeof(int);
    path_dest->points = malloc(len);
    memcpy(path_dest->points, path_src->points, len);
    path_dest->dlen = path_src->dlen;
    path_dest->plen = path_src->plen;
}

void free_paths(path_t *path, int n) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < path->dlen; j++) {
            free(&path[i].points[j]);
        }
    }
    free(path);
}

void free_points(pt_t *point, int n) {
    for (int i = 0; i < n; i++) {
        free(&point[i].data);
    }
    free(point);
}

void free_segmented(seg_t *seg) {
    free_paths(seg->paths, seg->n_seg);
    free_points(seg->points, seg->rsize);
    free(seg);
}

double get_nearest_score(double average_dist, double path_len) {
    return average_dist + path_len;
}

// void lines_pv2c(const int* lines, int lsize, int **cl) {

// }

// Converts PyVista line representation to a segmented line list in C
// Warning: ensure that the `lines` variable is null-terminated!
void seg_pv2c(double *pts, const int *lines, int lsize, seg_t *seg) {
    dbg_printf("calling seg_pv2c\n");
    
    int size = lines[0];
    int i = 0;
    seg->paths = NULL;
    seg->points = NULL;
    seg->n_seg = 0;
    seg->rsize = 0;

    print_seg(seg);

    dbg_printf("First line size: %d\n", size);

    dbg_printf("counting the number of segments\n");

    for (int i = 0; i < lsize; i++) {
        dbg_printf("%d ", lines[i]);
    }
    dbg_printf("\n");

    // count the number of segments
    while (size != 0) {
        seg->n_seg++;
        seg->rsize += size;
        i += size + 1;
        size = lines[i];
        dbg_printf("%d\n", size);
    }

    seg->rsize -= seg->n_seg;

    print_seg(seg);

    path_t *paths = (path_t *)calloc(seg->n_seg, sizeof(path_t));

    seg->points = (pt_t *)malloc(seg->rsize * sizeof(pt_t));

    size = lines[0];
    i = 0;
    int n = 0;

    dbg_printf("build new array and count path lengths\n");

    // build new segmented list and count path lengths
    while (size != 0) {
        paths[n].points = (int *)calloc(size, sizeof(int));
        paths[n].plen = 0;
        paths[n].dlen = size;

        dbg_printf("%d, %d\n", i, i+size);

        for (int j = i+1; j < i + size + 1; j++) {
            seg->points[lines[j]].seg_i = n;
            seg->points[lines[j]].data = (pts + lines[j] * DIM);

            paths[n].points[j - (i+1)] = lines[j];

            dbg_printf("%d %d, ", lines[j], pts[lines[j]]);

            paths[n].plen += 1;
            // j - i != 1 ? vec_dist(
            //     seg->points[lines[j]].data,
            //     seg->points[lines[j-1]].data
            // ) : 0; // actually compute distance
            // dbg_printf("seg_i: %d\n", seg->points[lines[j]].seg_i);
        }
        dbg_printf("\n");

        // for (int x = 0; x < seg->rsize; x++) {
        //     dbg_printf("%d ", seg->points[x].data);
        // }
        // dbg_printf("\n");

        dbg_printf("plen: %d\n", paths[n].plen);

        n++;
        i += size + 1;
        size = lines[i];
    }

    seg->paths = paths;

    dbg_printf("plen: %d\n", seg->paths[0].plen);

    print_seg(seg);

    dbg_printf("%ld\n", seg);

    dbg_printf("n_seg: %d\n", seg->n_seg);
}

int seg_c2pv(const seg_t *seg, int **lines) {
    print_seg(seg);

    dbg_printf("%d\n", (seg->rsize + seg->n_seg + 1));

    int asize = seg->rsize + seg->n_seg + 1;

    *lines = malloc(asize * sizeof(int));

    int *lp = *lines;

    dbg_printf("%d\n", *lines);

    dbg_printf("allocated space\n");

    int off = 0;

    for (int i = 0; i < seg->n_seg; i++) {
        lp[off] = seg->paths[i].dlen;
        dbg_printf("%d ", off);
        for (int j = 0; j < seg->paths[i].dlen; j++) {
            dbg_printf("%d ", seg->paths[i].points[j]);
            lp[off+j+1] = seg->paths[i].points[j];
        }
        off += seg->paths[i].dlen + 1;
        dbg_printf("\n");
    }

    print_arr(lp, (asize));

    // null terminate
    lp[asize - 1] = 0;

    return asize;
}

void skeletonize(const mesh_t *mesh, int *skeleton) {
    skeleton = (int *)calloc(mesh->vsize + 1, sizeof(int));

    // skeletonize mesh
}

void segment(const mesh_t *mesh, const int *lines, seg_t *seg) {
    seg->n_seg = 0;

    // segment mesh
}

void nearest_neighbors(const seg_t* ref, seg_t *tgt) {
    // paths on the target skeleton which already have a correspondence
    print_arr(ref->paths[0].points, ref->paths[0].dlen);
    print_arr(tgt->paths[0].points, tgt->paths[0].dlen);

    int *processed = (int *)calloc(tgt->n_seg, sizeof(int));

    int psize = 0;

    dbg_printf("TRAVERSE\n");

    for (int i = 0; i < ref->n_seg; i++) {
        double min_score = __DBL_MAX__;
        int min_path = 0;

        for (int j = 0; j < tgt->n_seg; j++) {
            if (processed[j])
                continue;
            
            double tdist = 0;
            double dist = 0;

            path_t r_path = ref->paths[i];

            path_t t_path = tgt->paths[j];

            printf("t_path plen: %f\n", t_path.plen);

            dbg_printf("FIND CLOSEST\n");
            
            // for each point in the reference path, find the nearest neighbor
            // in the opposing mesh and increment the total distance
            for (int rk = 0; rk < r_path.dlen; rk++) {
                // for (int x = 0; x < tgt->rsize; x++) {
                //     dbg_printf("%d ", tgt->points[x].data);
                // }
                // dbg_printf("\n");

                // dbg_printf("%d\n", r_path.points[rk]);
                // dbg_printf("%d\n", tgt->points[r_path.points[rk]]);
                // dbg_printf("%f\n", tgt->points[r_path.points[rk]].data[1]);

                for (int tk = 0; tk < t_path.dlen; tk++) {
                    // dbg_printf("call BLAS\n");
                    double d = vec_dist(tgt->points[t_path.points[tk]].data, 
                                        ref->points[r_path.points[rk]].data);

                    // dbg_printf("%f\n", d);
                    
                    if (d < dist)
                        dist = d;
                }

                dbg_printf("smallest distance: %f\n", dist);

                tdist += dist;
            }

            // score the target path's correspondence to the reference
            double score = get_nearest_score(
                tdist / ((double) r_path.dlen),
                fabs(t_path.plen - r_path.plen)
            );

            dbg_printf("plens: %f, %f\n", t_path.plen, r_path.plen);
            dbg_printf("tdist: %f\n", tdist);
            dbg_printf("distance: %f\n", tdist / ((double) r_path.dlen));
            dbg_printf("path length: %f\n", fabs(t_path.plen - r_path.plen));
            dbg_printf("score: %f\n", score);

            if (score < min_score) {
                min_score = score;
                min_path = j;
            }   
        }

        dbg_printf("min path: %d\n", min_path);

        processed[min_path] = i+1;
        psize += tgt->paths[min_path].dlen;
    }

    dbg_printf("SET PATHS\n");

    seg_t* new_seg = (seg_t *)calloc(1, sizeof(seg_t));

    // allocate path space for the number of paths in the reference
    new_seg->paths = (path_t *)malloc(ref->n_seg * sizeof(path_t));

    dbg_printf("n_seg: %d\n", ref->n_seg);

    new_seg->n_seg = ref->n_seg;

    int path_ind = 0;

    // add the maximum paths to the new seglist
    for (int i = 0; i < tgt->n_seg; i++) {
        dbg_printf("%d: %d\n", i, processed[i]);
        if (!processed[i])
            continue;

        path_cpy(&new_seg->paths[path_ind], &tgt->paths[i]);

        // print_arr(tgt->paths[path_ind].points, tgt->paths[path_ind].dlen);

        path_ind++;
    }

    // dbg_printf("n_seg: %d\n", new_seg->n_seg);
    // print_arr(new_seg->paths[0].points, new_seg->paths[0].dlen);

    // dbg_printf("points: %d\n", new_seg->points);

    dbg_printf("SET POINTS\n");

    print_seg(new_seg);

    // allocate space for only the total number of points in all maximum paths
    new_seg->points = (pt_t *)malloc(psize * sizeof(pt_t));

    new_seg->rsize = psize;

    dbg_printf("memory allocated\n");

    int point_ind = 0;

    // update points in target and add those with corresponding maximum paths
    for (int i = 0; i < ((int)tgt->rsize); i++) {
        int proc_val = processed[tgt->points[i].seg_i];
        if (!proc_val)
            continue;

        pt_cpy(&new_seg->points[point_ind], &tgt->points[i]);
        tgt->points[i].seg_i = proc_val-1;
    }

    dbg_printf("finished processing, freeing variables\n");

    // free stale memory ! THIS IS AN ISSUE
    // free_paths(tgt->paths, tgt->n_seg);
    // free_points(tgt->points, tgt->rsize);
    // free(tgt);

    dbg_printf("data freed\n");

    tgt->paths = new_seg->paths;
    tgt->points = new_seg->points;
    tgt->n_seg = new_seg->n_seg;
    tgt->rsize = new_seg->rsize;

    print_seg(tgt);

    dbg_printf("data set, exiting\n");
}

void c_find_skeleton(const mesh_t *mesh, seg_t *seg) {
    int *skeleton = NULL;

    skeletonize(mesh, skeleton);

    segment(mesh, skeleton, seg);
}

void c_nearest_neighbors(const seg_t* ref, const mesh_t *mesh, seg_t *tseg) {
    int *t_skeleton = NULL;

    skeletonize(mesh, t_skeleton);

    segment(mesh, t_skeleton, tseg);

    nearest_neighbors(ref, tseg);
}

void c_neighbors_from_segmented(const seg_t *ref, double *pts,
                                const int *skeleton, int ssize, seg_t *tseg) {
    seg_pv2c(pts, skeleton, ssize, tseg);

    nearest_neighbors(ref, tseg);
}

void c_seg_pv2c(double *pts, const int *lines, int lsize, seg_t *seg) {
    dbg_printf("[INFO]: converting to c skeleton\n");

    seg_pv2c(pts, lines, lsize, seg);

    dbg_printf("%ld\n", seg);
}

int main() {
    return 0;
}