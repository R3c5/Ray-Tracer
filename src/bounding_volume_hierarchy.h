#pragma once
#include "common.h"
#include <array>
#include <framework/ray.h>
#include <vector>
#include <queue>
#include <functional>

// Forward declaration.
struct Scene;

class BoundingVolumeHierarchy {
public:
    
    // Constructor. Receives the scene and builds the bounding volume hierarchy.
    BoundingVolumeHierarchy(Scene* pScene, const Features& features);

    // Return how many levels there are in the tree that you have constructed.
    [[nodiscard]] int numLevels() const;

    // Return how many leaf nodes there are in the tree that you have constructed.
    [[nodiscard]] int numLeaves() const;

    // Visual Debug 1: Draw the bounding boxes of the nodes at the selected level.
    void debugDrawLevel(int level);

    // Visual Debug 2: Draw the triangles of the i-th leaf
    void debugDrawLeaf(int leafIdx);

    // Return true if something is hit, returns false otherwise.
    // Only find hits if they are closer than t stored in the ray and the intersection
    // is on the correct side of the origin (the new t >= 0).
    bool intersect(Ray& ray, HitInfo& hitInfo, const Features& features) const;


    struct Node {
        //Type 0 = internal node;
        //Type 1 = leaf node
        bool type;
        int level = -1;
        AxisAlignedBox aabb;

        //For type 0 - contents will only contain 2 int representing the indices of the two child nodes
        //For type 1 - contents will contain the indices of both meshes and triangles within meshes in the following way
        //Meshes will be represented with negative numbers (and shifted by one so that 0 is reserved for triangles)
        //    For example: -1 will be mesh with index 0
        //Triangles will be represented with their normal index
        //    For example: 1 will be triangle with index 1
        std::vector<int> contents;
        
        //These are only used for SAH debug puproses
        int actual_split = -1;
        std::vector<std::pair<int, std::pair<int, int>>> debug_splits;

    };

    void createBVHTree(int parent, int level, const Features& features);
    
    //This is for debugging SAH
    void debugSAH(int level, int split);

    private:

    const int MAX_LEVELS = 100;
    int m_numLevels;
    int m_numLeaves;
    Scene* m_pScene;
    std::vector<Node> tree;


    bool basicTraversal(Ray& ray, HitInfo& hitInfo, int nodeInd, Vertex& v0Final, Vertex& v1Final, Vertex& v2Final, Features features, AxisAlignedBox& lastBox) const;
    
    #define pfn std ::pair<float, Node>
    bool traversalPQ(Ray& ray, HitInfo& hitInfo,
        std::priority_queue<pfn, std::vector<pfn>, std::function<bool(pfn, pfn)>>& pq,
        Vertex& v0Final, Vertex& v1Final, Vertex& v2Final, Features const& features) const;

    bool traversalVec(Ray& ray, HitInfo& hitInfo,
        std::vector<pfn>& pq,
        Vertex& v0Final, Vertex& v1Final, Vertex& v2Final, Features const& features) const;
};