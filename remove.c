#include <stdlib.h>

#include <cblas.h>

#define DIM 3

struct path {
    int *points;
    int plen;
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

double get_nearest_score(double average_dist, long path_len) {
    return average_dist + path_len;
}

void lines_pv2c(const int* lines, int lsize, int **cl) {

}

// Converts PyVista line representation to a segmented line list in C
void seg_pv2c(double **pts, const int *lines, int lsize, seg_t *seg) {
    int size = lines[0];
    int i = 0;
    seg->n_seg = 0;

    while (size != 0) {
        i += size;
        size = lines[i];
        seg->n_seg++;
    }

    path_t *paths = (path_t *)calloc(seg->n_seg, sizeof(path_t*));

    size = lines[0];
    i = 0;
    int n = 0;

    while (size != 0) {
        i += size;
        size = lines[i];
        paths[n].points = (int *)calloc(size, sizeof(int));
        paths[n].plen = 0;

        for (int j = i; j < i + size; j++) {
            seg->points[lines[j]].seg_i = n;
            seg->points[lines[j]].data = pts[lines[j]];

            paths[n].points[j - i] = lines[j];
            paths[n].plen += 1; // actually compute distance
        }

        n++;
    }

    seg->paths = paths;
}

void seg_c2pv(const seg_t *seg, int *lines) {
    for (int i = 0; i < seg->n_seg; i++) {
        lines[i*seg->n_seg] = seg->n_seg;
        for (int j = 0; j < seg->paths[i].dlen; j++) {
            lines[i*(seg->n_seg+1)+j] = seg->paths[i].points[j];
        }
    }
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
    int *processed = (int *)calloc(tgt->n_seg, sizeof(int));

    int psize = 0;

    for (int i = 0; i < ref->n_seg; i++) {
        double max_score = 0;
        int max_path = 0;

        for (int j = 0; j < tgt->n_seg; j++) {
            if (processed[j])
                continue;
            
            double tdist = 0;
            double dist = 0;

            path_t r_path = ref->paths[i];

            path_t t_path = tgt->paths[j];
            
            // for each point in the reference path, find the nearest neighbor
            // in the opposing mesh and increment the total distance
            for (int rk = 0; rk < r_path.dlen; rk++) {
                for (int tk = 0; tk < t_path.dlen; tk++) {
                    double d = vec_dist(tgt->points[t_path.points[tk]].data, 
                                        ref->points[r_path.points[rk]].data);
                    
                    if (d < dist)
                        dist = d;
                }

                tdist += dist;
            }

            // score the target path's correspondence to the reference
            double score = get_nearest_score(
                dist / ((double) r_path.dlen),
                t_path.plen
            );

            if (score > max_score) {
                max_score = score;
                max_path = j;
            }   
        }

        processed[max_path] = i;
        psize += tgt->paths[max_path].dlen;
    }

    seg_t* new_seg = (seg_t *)malloc(sizeof(seg_t));

    int path_ind = 0;

    // allocate path space for the number of paths in the reference
    new_seg->paths = (path_t *)malloc(ref->n_seg * sizeof(path_t));

    new_seg->n_seg = ref->n_seg;

    // add the maximum paths to the new seglist
    for (int i = 0; i < tgt->n_seg; i++) {
        if (!processed[i])
            continue;
        
        new_seg->paths[path_ind] = tgt->paths[i];
        path_ind++;
    }

    // allocate space for only the total number of points in all maximum paths
    new_seg->points = (pt_t *)malloc(psize * sizeof(pt_t));

    int point_ind = 0;

    // update points in target and add those with corresponding maximum paths
    for (int i = 0; i < tgt->rsize; i++) {
        int proc_val = processed[tgt->points[i].seg_i];
        if (!proc_val)
            continue;
            
        tgt->points[i].seg_i = proc_val;
        new_seg->points[point_ind] = tgt->points[i];
    }

    // free stale memory
    free(tgt->paths);
    free(tgt->points);
    free(tgt);

    tgt = new_seg;
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

void c_neighbors_from_segmented(const seg_t *ref, double **pts,
                                const int *skeleton, int ssize, seg_t *tseg) {
    seg_pv2c(pts, skeleton, ssize, tseg);

    nearest_neighbors(ref, tseg);
}

void c_seg_pv2c(double **pts, const int *lines, int lsize, seg_t *seg) {
    seg_pv2c(pts, lines, lsize, seg);
}

int main() {
    return 0;
}