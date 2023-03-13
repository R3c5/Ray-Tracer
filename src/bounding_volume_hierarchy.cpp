#include "bounding_volume_hierarchy.h"
#include "draw.h"
#include "intersect.h"
#include "scene.h"
#include "texture.h"
#include "interpolate.h"
#include <glm/glm.hpp>
#include <iostream>
#include <chrono>
#include <functional>


//Choose from the following options by commenting lines that are not wanted
//and uncommenting the ones that are prefered.

//Timer for creation of BVH
//#define creation_time

// Traversal type
//#define basic_traversal
//#define traversal_pq
#define traversal_vec

//Timer for traversal for each individual ray
//#define traversal_debug_time

//Number of bins used for SAH
#define numBins 5
//Enable Visual Debug for SAH
#define SAHdebug


AxisAlignedBox calculateAABB(Scene* pscene, std::vector<int> const& contents)
{

    glm::vec3 minB = glm::vec3(std::numeric_limits<float>::max());
    glm::vec3 maxB = -glm::vec3(std::numeric_limits<float>::max());
    //First element in vector should be a mesh index
    int meshIndex = -1;
    for (int i : contents) {
        if (i < 0) {
            meshIndex = - i - 1; 
            continue;
        }
        glm::uvec3& tri = pscene->meshes[meshIndex].triangles[i];
        glm::vec3& vert1 = pscene->meshes[meshIndex].vertices[tri.x].position;
        glm::vec3& vert2 = pscene->meshes[meshIndex].vertices[tri.y].position;
        glm::vec3& vert3 = pscene->meshes[meshIndex].vertices[tri.z].position;

        float maxBX = std::max(std::max(maxB.x, vert1.x), std::max(vert2.x, vert3.x));
        float maxBY = std::max(std::max(maxB.y, vert1.y), std::max(vert2.y, vert3.y));
        float maxBZ = std::max(std::max(maxB.z, vert1.z), std::max(vert2.z, vert3.z));
        float minBX = std::min(std::min(minB.x, vert1.x), std::min(vert2.x, vert3.x));
        float minBY = std::min(std::min(minB.y, vert1.y), std::min(vert2.y, vert3.y));
        float minBZ = std::min(std::min(minB.z, vert1.z), std::min(vert2.z, vert3.z));

        minB = glm::vec3(minBX, minBY, minBZ);
        maxB = glm::vec3(maxBX, maxBY, maxBZ);
    }
    
    return AxisAlignedBox { minB, maxB };
}
AxisAlignedBox calculateAABBFromPairs(Scene* pscene, std::vector<std::pair<int,int>> const& contents)
{

    glm::vec3 minB = glm::vec3(std::numeric_limits<float>::max());
    glm::vec3 maxB = -glm::vec3(std::numeric_limits<float>::max());
    // First element in vector should be a mesh index
    int meshIndex = -1;
    for (std::pair<int,int> p : contents) {
        
        glm::uvec3& tri = pscene->meshes[p.first].triangles[p.second];
        glm::vec3& vert1 = pscene->meshes[p.first].vertices[tri.x].position;
        glm::vec3& vert2 = pscene->meshes[p.first].vertices[tri.y].position;
        glm::vec3& vert3 = pscene->meshes[p.first].vertices[tri.z].position;

        float maxBX = std::max(std::max(maxB.x, vert1.x), std::max(vert2.x, vert3.x));
        float maxBY = std::max(std::max(maxB.y, vert1.y), std::max(vert2.y, vert3.y));
        float maxBZ = std::max(std::max(maxB.z, vert1.z), std::max(vert2.z, vert3.z));
        float minBX = std::min(std::min(minB.x, vert1.x), std::min(vert2.x, vert3.x));
        float minBY = std::min(std::min(minB.y, vert1.y), std::min(vert2.y, vert3.y));
        float minBZ = std::min(std::min(minB.z, vert1.z), std::min(vert2.z, vert3.z));

        minB = glm::vec3(minBX, minBY, minBZ);
        maxB = glm::vec3(maxBX, maxBY, maxBZ);
    }

    return AxisAlignedBox { minB, maxB };
}

glm::vec3 calculateCentroidOfTriangle(glm::vec3 &v1, glm::vec3 &v2, glm::vec3 &v3) {
    float x = (v1.x + v2.x + v3.x) / 3.0f;
    float y = (v1.y + v2.y + v3.y) / 3.0f;
    float z = (v1.z + v2.z + v3.z) / 3.0f;

    return glm::vec3 (x, y, z);
}
glm::vec3 getCentroidFromIndex(Mesh &mesh, int triInd)
{

    glm::uvec3& tri = mesh.triangles[triInd];
    glm::vec3& vert1 = mesh.vertices[tri.x].position;
    glm::vec3& vert2 = mesh.vertices[tri.y].position;
    glm::vec3& vert3 = mesh.vertices[tri.z].position;

    return calculateCentroidOfTriangle(vert1, vert2, vert3);
}

void splitByMedianTriangle(Scene* pscene,
    std::vector<int> const& contents,
    std::vector<int> &left,
    std::vector<int> &right,
    int const level)
{
    int axis = (level - 1) % 3;
    glm::vec3 axisCoord = (axis == 0) ? glm::vec3(1, 0, 0) : 
                          (axis == 1) ? glm::vec3(0, 1, 0) : glm::vec3(0, 0, 1);

    //vector contains the coordinate of the centroid of the triangle on the axis, and a pair of meshes index and triangle index
    std::vector<std::pair<float, std::pair<int, int>>> triVec;
    
    int meshIndex = -1;
     for (int i : contents) {
        if (i < 0) {
            meshIndex = -i - 1;
            continue;
        }
        glm::vec3 triCentroid = getCentroidFromIndex(pscene->meshes[meshIndex], i);
        float centroid = glm::dot(triCentroid, axisCoord);
        triVec.push_back(std::pair(centroid,std::pair(meshIndex, i)));
    }

    //Calculate the median triangle
    //Syntax found on https://en.cppreference.com/w/cpp/algorithm/nth_element
    auto median = triVec.begin() + triVec.size() / 2;
    std::nth_element(triVec.begin(), median, triVec.end());

    // Do the split
    left = {};
    right = {};
    //Split the elements to the left child
    int currentMeshIndex = 0;
    int i;
    for (i = 0; i < (triVec.size() + 1) / 2; i++) {
        int meshInd = - triVec[i].second.first - 1;
        if (meshInd != currentMeshIndex) {
            left.push_back(meshInd);
            currentMeshIndex = meshInd;
        }
        left.push_back(triVec[i].second.second);
    }
    
    //Split the elements to the right child
    currentMeshIndex = 0;
    for (i = (triVec.size() + 1) / 2; i < triVec.size(); i++) {
        int meshInd = -triVec[i].second.first - 1;
        if (meshInd != currentMeshIndex) {
            right.push_back(meshInd);
            currentMeshIndex = meshInd;
        }
        right.push_back(triVec[i].second.second);
    }
}

float calculateVolume(AxisAlignedBox aabb)
{
    return abs(aabb.upper.x - aabb.lower.x) * abs(aabb.upper.y - aabb.lower.y) * abs(aabb.upper.z - aabb.lower.z);
}

void distributeSplit(std::vector<int>& vec, std::vector<std::pair<int, int>> & split)
{
    vec = {};
    int prevIndex = 0, curInd;
    for (std::pair<int, int> p : split) {
        curInd = -p.first - 1;
        if (curInd != prevIndex) {
            prevIndex = curInd;
            vec.push_back(prevIndex);
        }
        vec.push_back(p.second);
    }
}

void splitBySAH(Scene* pscene,
    BoundingVolumeHierarchy::Node &parent,
    std::vector<int>& left,
    std::vector<int>& right) {

    float sdx = abs(parent.aabb.upper.x - parent.aabb.lower.x);
    float sdy = abs(parent.aabb.upper.y - parent.aabb.lower.y);
    float sdz = abs(parent.aabb.upper.z - parent.aabb.lower.z);
    int axis;
    float len, start;
    if (sdx >= sdy && sdx >= sdz) {
        axis = 0;
        len = sdx;
        start = std::min(parent.aabb.lower.x, parent.aabb.upper.x);
    } else if (sdy >= sdx && sdy >= sdz ) {
        axis = 1;
        len = sdy;
        start = std::min(parent.aabb.lower.y, parent.aabb.upper.y);
    } else {
        axis = 2;
        len = sdz;
        start = std::min(parent.aabb.lower.z, parent.aabb.upper.z);
    }

    glm::vec3 axisCoord = (axis == 0) ? glm::vec3(1, 0, 0) : (axis == 1) ? glm::vec3(0, 1, 0)
                                                                         : glm::vec3(0, 0, 1);
    

    std::vector<std::pair<int, std::pair<int, int>>> bins;
    
    int meshIndex = -1;
    for (int i : parent.contents) {
        if (i < 0) {
            meshIndex = -i - 1;
            continue;
        }
        glm::vec3 triCentroid = getCentroidFromIndex(pscene->meshes[meshIndex], i);
        float centroid = glm::dot(triCentroid, axisCoord);

        int bin = (centroid - start) * numBins / len;
        if (bin == numBins)
            bin--;
        bins.push_back(std::pair(bin, std::pair(meshIndex, i)));
    }
    std::sort(bins.begin(), bins.end());

    #ifdef SAHdebug
        parent.debug_splits = bins;
    #endif

    left = {};
    right = {};

    std::vector<std::pair<int, int>> lSplit;
    std::vector<std::pair<int, int>> rSplit;
 
    float cost_min = bins.size() * calculateVolume(parent.aabb); // std::numeric_limits<float>::max();
    //int pos = 0;
    int min_size_dif = bins.size();
    for (int split = 0; split < numBins -1; split++) {
        lSplit = {};
        rSplit = {};
        for (int i = 0; i < bins.size(); i++) {
            if (bins[i].first <= split) {
                lSplit.push_back(bins[i].second);
            } else {
                rSplit.push_back(bins[i].second);
            }
        }
        if (lSplit.size() == 0 || rSplit.size() == 0) {
            continue;
        }
        float cost = lSplit.size() * calculateVolume(calculateAABBFromPairs(pscene, lSplit)) + rSplit.size() * calculateVolume(calculateAABBFromPairs(pscene, rSplit));
        
        if (abs(cost - cost_min) < 0.00001) {
            int dif = lSplit.size() - rSplit.size();
            int cur_dif = std::abs(dif);
            if (min_size_dif > cur_dif) {
                min_size_dif = cur_dif;
                #ifdef SAHdebug
                    parent.actual_split = split;
                #endif

                distributeSplit(left, lSplit);
                distributeSplit(right, rSplit);
            }
        }
        
        if (cost < cost_min) {
            cost_min = cost;

            #ifdef SAHdebug
                parent.actual_split = split;
            #endif

            distributeSplit(left, lSplit);
            distributeSplit(right, rSplit);
        }

    }

        

}

void BoundingVolumeHierarchy::createBVHTree(int parent, int level, const Features& features)
{
   
    this->m_numLevels = std::max(numLevels(), level);

    if (this->tree[parent].type == 0)
        return;
    
    if (level >= MAX_LEVELS || this->tree[parent].contents.size() <= 2 /*checks if contents only has one triangle (index of mesh and triangle)*/) {
        this->m_numLeaves++;
        return;
    }


    std::vector<int> leftCont = {};
    std::vector<int> rightCont = {};
    if (!features.extra.enableBvhSahBinning) {
        splitByMedianTriangle(m_pScene, this->tree[parent].contents, leftCont, rightCont, level);
    } else {
        splitBySAH(m_pScene, this->tree[parent], leftCont, rightCont);
    }
    if (leftCont.size() == 0 || rightCont.size() == 0)
        return;
    //Create left child
    int leftIndex = this->tree.size();
    Node leftNode;
    leftNode.type = 1;
    leftNode.level = level;
    leftNode.contents = leftCont;
    leftNode.aabb = calculateAABB(m_pScene, leftCont);
    this->tree.push_back(leftNode);

    //Create right child.
    int rightIndex = leftIndex + 1;
    Node rightNode;
    rightNode.type = 1;
    rightNode.level = level;
    rightNode.contents = rightCont;
    rightNode.aabb = calculateAABB(m_pScene, rightCont);
    this->tree.push_back(rightNode);

    //Update parent node
    this->tree[parent].type = 0;
    this->tree[parent].contents = { leftIndex, rightIndex };

    this->createBVHTree(leftIndex, level + 1, features);
    this->createBVHTree(rightIndex, level + 1, features);

}

BoundingVolumeHierarchy::BoundingVolumeHierarchy(Scene* pScene, const Features& features)
    : m_pScene(pScene)
{

    #ifdef creation_time
    auto t1 = std::chrono::high_resolution_clock::now();
    #endif
    
    this->m_pScene = pScene;
    this->m_numLevels = 0;
    this->m_numLeaves = 0;
    // TODO: implement BVH construction.
    std::vector<int> allIndices;
    for (int i = 0; i < pScene->meshes.size(); i++) {
        allIndices.push_back(-i - 1);
        Mesh mesh = pScene->meshes[i];
    for (int j = 0; j < mesh.triangles.size(); j++) {
            allIndices.push_back(j);
        }
    }
     
    Node root;
    root.type = 1; //leaf for now
    root.level = 0;
    root.contents = allIndices;
    root.aabb = calculateAABB(pScene, allIndices);

    this->tree.push_back(root);

    //Start with root at index 0 and create children on level 1
    createBVHTree(0, 1, features);

    #ifdef creation_time
    auto t2 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> delta = t2 - t1;
    std::cout << "Time elapsed: " << delta.count() << " ms\n";
    #endif
}


// Return the depth of the tree that you constructed. This is used to tell the
// slider in the UI how many steps it should display for Visual Debug 1.
int BoundingVolumeHierarchy::numLevels() const
{
    return m_numLevels;
}

// Return the number of leaf nodes in the tree that you constructed. This is used to tell the
// slider in the UI how many steps it should display for Visual Debug 2.
int BoundingVolumeHierarchy::numLeaves() const
{
    return m_numLeaves;
}

// Use this function to visualize your BVH. This is useful for debugging. Use the functions in
// draw.h to draw the various shapes. We have extended the AABB draw functions to support wireframe
// mode, arbitrary colors and transparency.
void BoundingVolumeHierarchy::debugDrawLevel(int level)
{
    // Draw the AABB as a transparent green box.
    //AxisAlignedBox aabb{ glm::vec3(-0.05f), glm::vec3(0.05f, 1.05f, 1.05f) };
    //drawShape(aabb, DrawMode::Filled, glm::vec3(0.0f, 1.0f, 0.0f), 0.2f);

    // Draw the AABB as a (white) wireframe box.
    // AxisAlignedBox aabb { glm::vec3(0.0f), glm::vec3(0.0f, 1.05f, 1.05f) };
    //drawAABB(aabb, DrawMode::Wireframe);
    // drawAABB(aabb, DrawMode::Filled, glm::vec3(0.05f, 1.0f, 0.05f), 0.1f);

    for (int i = 0; i < tree.size(); i++) {
        Node node = tree[i];
        if (node.level != level)
            continue;
        drawAABB(node.aabb, DrawMode::Filled, glm::vec3(0.05f, 1.0f, 0.05f), 0.1f);
    }
}



// Use this function to visualize your leaf nodes. This is useful for debugging. The function
// receives the leaf node to be draw (think of the ith leaf node). Draw the AABB of the leaf node and all contained triangles.
// You can draw the triangles with different colors. NoteL leafIdx is not the index in the node vector, it is the
// i-th leaf node in the vector.
void BoundingVolumeHierarchy::debugDrawLeaf(int leafIdx)
{
    // Draw the AABB as a transparent green box.
    //AxisAlignedBox aabb{ glm::vec3(-0.05f), glm::vec3(0.05f, 1.05f, 1.05f) };
    //drawShape(aabb, DrawMode::Filled, glm::vec3(0.0f, 1.0f, 0.0f), 0.2f);

    // Draw the AABB as a (white) wireframe box.
    // AxisAlignedBox aabb { glm::vec3(0.0f), glm::vec3(0.0f, 1.05f, 1.05f) };
    //drawAABB(aabb, DrawMode::Wireframe);
    // drawAABB(aabb, DrawMode::Filled, glm::vec3(0.05f, 1.0f, 0.05f), 0.1f);

    int leafNo = 0;
    for (int ind = 0; ind < tree.size(); ind++) {
        if (tree[ind].type == 1) {
            leafNo++;

            if (leafNo == leafIdx) {
                Node node = tree[ind];
                drawAABB(node.aabb, DrawMode::Wireframe, glm::vec3(0.05f, 0.05f, 1.0f), 0.5f);

                int meshIndex = -1;
                for (int i : node.contents) {
                    if (i < 0) {
                        meshIndex = -i - 1;
                        continue;
                    }
                    
                    Mesh& mesh = m_pScene->meshes[meshIndex];
                    glm::uvec3& tri = mesh.triangles[i];
                    Vertex vert1 = mesh.vertices[tri.x];
                    Vertex vert2 = mesh.vertices[tri.y];
                    Vertex vert3 = mesh.vertices[tri.z];
                    drawTriangle(vert1, vert2, vert3);
                }
                return;
            }
        }
    }
    // once you find the leaf node, you can use the function drawTriangle (from draw.h) to draw the contained primitives
}


void BoundingVolumeHierarchy::debugSAH(int level, int split)
{
    
    for (int i = 0; i < this->tree.size(); i++) {
        Node node = this->tree[i];
        if (node.level != level)
            continue;

        drawAABB(node.aabb, DrawMode::Wireframe, glm::vec3(1.0f, 1.0f, 1.0f), 1.0f);
        

        #ifdef SAHdebug

        if (node.actual_split == -1)
            continue;

        std::vector<std::pair<int, int>> lSplit;
        std::vector<std::pair<int, int>> rSplit;
        for (int i = 0; i < node.debug_splits.size(); i++) {
            if (node.debug_splits[i].first <= split) {
                lSplit.push_back(node.debug_splits[i].second);
            } else {
                rSplit.push_back(node.debug_splits[i].second);
            }
        }
        glm::vec3 color = glm::vec3(1.0f, 0.05f, 0.05f);
        if (split == node.actual_split)
            color = glm::vec3(0.05f, 1.0f, 0.05f);
        if (lSplit.size() > 0) {
            AxisAlignedBox box = calculateAABBFromPairs(m_pScene, lSplit);
            drawAABB(box, DrawMode::Filled, color, 0.5f);
        }
        if (rSplit.size() > 0) {
            AxisAlignedBox box = calculateAABBFromPairs(m_pScene, rSplit);
            drawAABB(box, DrawMode::Filled, color, 0.5f);
        }
        #endif

    }
}

bool BoundingVolumeHierarchy::basicTraversal(Ray& ray, HitInfo& hitInfo, int nodeInd, Vertex& v0Final, Vertex& v1Final, Vertex& v2Final, Features features, AxisAlignedBox& lastBox) const
{

    Node node = this->tree[nodeInd];
    float t = ray.t;
    bool hit = intersectRayWithShape(node.aabb, ray);
    ray.t = t;
    // If box is not hit return false
    if (!hit) {
        if (features.enableExtraBoxes)
        drawAABB(node.aabb, DrawMode::Wireframe, glm::vec3(1.0f, 0.05f, 0.05f), 0.2f);
        return false;
    }
       
    hit = false;
    // Leaf node - check all triangles

    if (node.type == 1) {
        int meshInd = -1;
        for (int i : node.contents) {
            if (i < 0) {
                meshInd = -i - 1;
                continue;
            }
            Mesh& mesh = this->m_pScene->meshes[meshInd];
            const auto& tri = mesh.triangles[i];
            const auto& v0 = mesh.vertices[tri.x];
            const auto& v1 = mesh.vertices[tri.y];
            const auto& v2 = mesh.vertices[tri.z];

            if (intersectRayWithTriangle(v0.position, v1.position, v2.position, ray, hitInfo)) {
                hitInfo.material = mesh.material;
                v0Final = v0;
                v1Final = v1;
                v2Final = v2;
                if (features.enableNormalInterp) {
                    glm::vec3 bary = computeBarycentricCoord(v0Final.position, v1Final.position, v2Final.position, ray.origin + ray.t * ray.direction);
                    glm::vec3 interpolatedNormal = interpolateNormal(v0Final.normal, v1Final.normal, v2Final.normal, bary);
                    hitInfo.normal = interpolatedNormal;
                    if (glm::dot(glm::normalize(ray.direction), interpolatedNormal) > 0)
                        hitInfo.normal = -1.0f * interpolatedNormal;
                } else
                    hitInfo.normal = glm::normalize(glm::cross(v1.position - v0.position, v2.position - v0.position));
                lastBox = node.aabb;
                hit = true;
            }
        }
        return hit;
    }
    drawAABB(node.aabb, DrawMode::Wireframe, glm::vec3(0.05f, 1.0f, 0.05f), 0.1f);
    int leftInd = node.contents[0];
    int rightInd = leftInd + 1;
    hit = basicTraversal(ray, hitInfo, leftInd, v0Final, v1Final, v2Final, features, lastBox);
    hit |= basicTraversal(ray, hitInfo, rightInd, v0Final, v1Final, v2Final, features, lastBox);
    return hit;
}
#define pfn std ::pair<float, Node>
#define comparator [](pfn a, pfn b){ return a.first < b.first; }
bool BoundingVolumeHierarchy::traversalPQ(Ray& ray, HitInfo& hitInfo,
    std::priority_queue<pfn, std::vector<pfn>, std::function<bool(pfn, pfn)>>& pq,
    Vertex& v0Final, Vertex& v1Final, Vertex& v2Final, Features const& features) const
{
    if (pq.empty()) {
        return false;
    }
    
    std::pair<float, Node> pair = pq.top();
    pq.pop();
    Node node = pair.second;
    
    bool hit = false;
    //Leaf node - check all triangles
    if (node.type == 1) {
        int meshInd = -1;
        for (int i : node.contents) {
            if (i < 0) {
                meshInd = -i - 1;
                continue;
            }
            Mesh& mesh = this->m_pScene->meshes[meshInd];
            const auto& tri = mesh.triangles[i];
            const auto& v0 = mesh.vertices[tri.x];
            const auto& v1 = mesh.vertices[tri.y];
            const auto& v2 = mesh.vertices[tri.z];
            if (intersectRayWithTriangle(v0.position, v1.position, v2.position, ray, hitInfo)) {
                hitInfo.material = mesh.material;
                v0Final = v0;
                v1Final = v1;
                v2Final = v2;
                if (features.enableNormalInterp) {
                    glm::vec3 bary = computeBarycentricCoord(v0Final.position, v1Final.position, v2Final.position, ray.origin + ray.t * ray.direction);
                    glm::vec3 interpolatedNormal = interpolateNormal(v0Final.normal, v1Final.normal, v2Final.normal, bary);
                    hitInfo.normal = interpolatedNormal;
                    if (glm::dot(glm::normalize(ray.direction), interpolatedNormal) > 0)
                        hitInfo.normal = -1.0f * interpolatedNormal;
                } else
                    hitInfo.normal = glm::normalize(glm::cross(v1.position - v0.position, v2.position - v0.position));
                hit = true;
            }
        }
        return hit | traversalPQ(ray, hitInfo, pq, v0Final, v1Final, v2Final, features);
    }

    //Interior node
    int leftInd = node.contents[0];
    int rightInd = leftInd + 1;
    float t = ray.t;
    bool lh = intersectRayWithShape(tree[leftInd].aabb, ray);
    float lt = ray.t;
    ray.t = t;
    bool rh = intersectRayWithShape(tree[rightInd].aabb, ray);
    float rt = ray.t;
    ray.t = t;
    
    if (lh) {
        pq.push(std::pair(lt, tree[leftInd]));
        drawAABB(tree[leftInd].aabb, DrawMode::Wireframe, glm::vec3(0.05f, 1.0f, 0.05f), 0.1f);
    } else {
        if (features.enableExtraBoxes)
        drawAABB(tree[leftInd].aabb, DrawMode::Wireframe, glm::vec3(1.0f, 0.05f, 0.05f), 0.1f);
    }
    if (rh) {
        pq.push(std::pair(rt, tree[rightInd]));
        drawAABB(tree[rightInd].aabb, DrawMode::Wireframe, glm::vec3(0.05f, 1.0f, 0.05f), 0.1f);
    } else {
        if (features.enableExtraBoxes)
        drawAABB(tree[rightInd].aabb, DrawMode::Wireframe, glm::vec3(1.0f, 0.05f, 0.05f), 0.1f);
    }
    return traversalPQ(ray, hitInfo, pq, v0Final, v1Final, v2Final,features);
}

#define comparatorDec [](pfn a, pfn b) { return a.first > b.first; }
bool BoundingVolumeHierarchy::traversalVec(Ray& ray, HitInfo& hitInfo,
    std::vector<pfn>& vec, Vertex& v0Final, Vertex& v1Final, Vertex& v2Final, Features const& features) const
{
    if (vec.empty()) {
        return false;
    }
   // std::nth_element(vec.begin(), vec.end() - 1, vec.end(), comparatorDec);
    std::pair<float, Node> pair = vec.back();
    vec.pop_back();
    Node node = pair.second;

    bool hit = false;
    // Leaf node - check all triangles
    if (node.type == 1) {
        int meshInd = -1;
        for (int i : node.contents) {
            if (i < 0) {
                meshInd = -i - 1;
                continue;
            }
            Mesh& mesh = this->m_pScene->meshes[meshInd];
            const auto& tri = mesh.triangles[i];
            const auto& v0 = mesh.vertices[tri.x];
            const auto& v1 = mesh.vertices[tri.y];
            const auto& v2 = mesh.vertices[tri.z];
            if (intersectRayWithTriangle(v0.position, v1.position, v2.position, ray, hitInfo)) {
                hitInfo.material = mesh.material;
                v0Final = v0;
                v1Final = v1;
                v2Final = v2;
                if (features.enableNormalInterp) {
                    glm::vec3 bary = computeBarycentricCoord(v0Final.position, v1Final.position, v2Final.position, ray.origin + ray.t * ray.direction);
                    glm::vec3 interpolatedNormal = interpolateNormal(v0Final.normal, v1Final.normal, v2Final.normal, bary);
                    hitInfo.normal = interpolatedNormal;
                    if (glm::dot(glm::normalize(ray.direction), interpolatedNormal) > 0)
                        hitInfo.normal = -1.0f * interpolatedNormal;
                } else
                    hitInfo.normal = glm::normalize(glm::cross(v1.position - v0.position, v2.position - v0.position));
                hit = true;
            }
        }
        return hit | traversalVec(ray, hitInfo, vec, v0Final, v1Final, v2Final, features);
    }

    // Interior node
    int leftInd = node.contents[0];
    int rightInd = leftInd + 1;
    float t = ray.t;
    bool lh = intersectRayWithShape(tree[leftInd].aabb, ray);
    float lt = ray.t;
    ray.t = t;
    bool rh = intersectRayWithShape(tree[rightInd].aabb, ray);
    float rt = ray.t;
    ray.t = t;

    if (lh) {
        vec.push_back(std::pair(lt, tree[leftInd]));
        drawAABB(tree[leftInd].aabb, DrawMode::Wireframe, glm::vec3(0.05f, 1.0f, 0.05f), 0.1f);
    } else {
        if (features.enableExtraBoxes)
        drawAABB(tree[leftInd].aabb, DrawMode::Wireframe, glm::vec3(1.0f, 0.05f, 0.05f), 0.1f);
    }
    if (rh) {
        vec.push_back(std::pair(rt, tree[rightInd]));
        drawAABB(tree[rightInd].aabb, DrawMode::Wireframe, glm::vec3(0.05f, 1.0f, 0.05f), 0.1f);
    } else {
        if (features.enableExtraBoxes)
        drawAABB(tree[rightInd].aabb, DrawMode::Wireframe, glm::vec3(1.0f, 0.05f, 0.05f), 0.1f);
    }
    std::sort(vec.begin(), vec.end(), comparatorDec);
    return traversalVec(ray, hitInfo, vec, v0Final, v1Final, v2Final, features);
}


// Return true if something is hit, returns false otherwise. Only find hits if they are closer than t stored
// in the ray and if the intersection is on the correct side of the origin (the new t >= 0). Replace the code
// by a bounding volume hierarchy acceleration structure as described in the assignment. You can change any
// file you like, including bounding_volume_hierarchy.h.
bool BoundingVolumeHierarchy::intersect(Ray& ray, HitInfo& hitInfo, const Features& features) const
{
    Vertex v0Final;
    Vertex v1Final;
    Vertex v2Final;
    // If BVH is not enabled, use the naive implementation.
    if (!features.enableAccelStructure) {
        bool hit = false;
        
        // Intersect with all triangles of all meshes.
        for (const auto& mesh : m_pScene->meshes) {
            for (const auto& tri : mesh.triangles) {
                const auto& v0 = mesh.vertices[tri[0]];
                const auto& v1 = mesh.vertices[tri[1]];
                const auto& v2 = mesh.vertices[tri[2]];
                if (intersectRayWithTriangle(v0.position, v1.position, v2.position, ray, hitInfo)) {
                    hitInfo.material = mesh.material;
                    v0Final = v0;
                    v1Final = v1;
                    v2Final = v2;
                    if (features.enableNormalInterp) {
                        glm::vec3 bary = computeBarycentricCoord(v0Final.position, v1Final.position, v2Final.position, ray.origin + ray.t * ray.direction);
                        glm::vec3 interpolatedNormal = interpolateNormal(v0Final.normal, v1Final.normal, v2Final.normal, bary);
                        hitInfo.normal = interpolatedNormal;
                        if (glm::dot(glm::normalize(ray.direction), interpolatedNormal) > 0)
                            hitInfo.normal = -1.0f * interpolatedNormal;
                    } else
                        hitInfo.normal = glm::normalize(glm::cross(v1.position - v0.position, v2.position - v0.position));
                    hit = true;
                }
            }
        }

        // Intersect with spheres.
        for (const auto& sphere : m_pScene->spheres)
            hit |= intersectRayWithShape(sphere, ray, hitInfo);

        if (features.enableTextureMapping) {
            glm::vec3 bary = computeBarycentricCoord(v0Final.position, v1Final.position, v2Final.position, ray.origin + ray.t * ray.direction);
            glm::vec2 texCoord = interpolateTexCoord(v0Final.texCoord, v1Final.texCoord, v2Final.texCoord, bary);
            hitInfo.texCoord = texCoord;
            if (hitInfo.material.kdTexture) {
                hitInfo.material.kd = acquireTexel(*hitInfo.material.kdTexture.get(), texCoord, features);
            }
        }

        if (features.enableNormalInterp == true) {
            glm::vec3 bary = computeBarycentricCoord(v0Final.position, v1Final.position, v2Final.position, ray.origin + ray.t * ray.direction);
            glm::vec3 interpolatedNormal = interpolateNormal(v0Final.normal, v1Final.normal, v2Final.normal, bary);
            Ray n0, n1, n2, nIterpolated;
            n0.origin = v0Final.position;
            n0.direction = v0Final.normal;
            n0.t = 1.0f;
            n1.origin = v1Final.position;
            n1.direction = v1Final.normal;
            n1.t = 1.0f;
            n2.origin = v2Final.position;
            n2.direction = v2Final.normal;
            n2.t = 1.0f;
            nIterpolated.origin = ray.origin + ray.t * ray.direction;
            nIterpolated.direction = interpolatedNormal;
            if (glm::dot(glm::normalize(ray.direction), interpolatedNormal) > 0)
                nIterpolated.direction = -1.0f * interpolatedNormal;
            nIterpolated.t = 1.0f;
            drawRay(n0);
            drawRay(n1);
            drawRay(n2);
            drawRay(nIterpolated);
        }

        return hit;
    } else {
        // TODO: implement here the bounding volume hierarchy traversal.
        // Please note that you should use `features.enableNormalInterp` and `features.enableTextureMapping`
        // to isolate the code that is only needed for the normal interpolation and texture mapping features.


        //Start timing
        #ifdef traversal_debug_time
        auto t1 = std::chrono::high_resolution_clock::now();
        #endif

        
        bool hit = false;

       // Basic traversal
        #ifdef basic_traversal
        AxisAlignedBox lastBox = tree[0].aabb;
        float prevT = ray.t;
        if (!intersectRayWithShape(tree[0].aabb, ray)) {
            if (features.enableExtraBoxes)
            drawAABB(tree[0].aabb, DrawMode::Wireframe, glm::vec3(1.0f, 0.05f, 0.05f), 0.2f);
            return false;
        }
        ray.t = prevT;
        hit = this->basicTraversal(ray, hitInfo, 0, v0Final, v1Final, v2Final, features, lastBox);
        drawAABB(lastBox, DrawMode::Wireframe, glm::vec3(0.05f, 1.0f, 0.05f), 0.2f);
        #endif

        //Traversal with Priority Queue
        #ifdef traversal_pq
            float prevT = ray.t;

            if (!intersectRayWithShape(tree[0].aabb, ray)) {
                if (features.enableExtraBoxes)
                drawAABB(tree[0].aabb, DrawMode::Wireframe, glm::vec3(1.0f, 0.05f, 0.05f), 0.2f);
                return false;
            }
            drawAABB(tree[0].aabb, DrawMode::Wireframe, glm::vec3(0.05f, 1.0f, 0.05f), 0.2f);
            std::priority_queue<pfn, std::vector<pfn>,
                std::function<bool(pfn, pfn)>>
                pq(comparator);
            pq.push(std::pair(ray.t, tree[0]));

            ray.t = prevT;
            hit = traversalPQ(ray, hitInfo, pq, v0Final, v1Final, v2Final, features);
        #endif

        // Traversal with vector
        #ifdef traversal_vec
            float prevT = ray.t;

            if (!intersectRayWithShape(tree[0].aabb, ray)) {
                if (features.enableExtraBoxes)
                drawAABB(tree[0].aabb, DrawMode::Wireframe, glm::vec3(1.0f, 0.05f, 0.05f), 0.2f);
                return false;
            }
            drawAABB(tree[0].aabb, DrawMode::Wireframe, glm::vec3(0.05f,1.0f, 0.05f), 0.2f);
            std::vector<pfn> vec;
            vec.push_back(std::pair(ray.t, tree[0]));

            ray.t = prevT;
            hit = traversalVec(ray, hitInfo, vec, v0Final, v1Final, v2Final, features);
        #endif
        if (hit)
            drawTriangleWithColor(v0Final.position, v1Final.position, v2Final.position, glm::vec3 { 0.0f, 0.0f, 1.0f }, 0.1f);

        // Intersect with spheres.
        for (const auto& sphere : m_pScene->spheres)
            hit |= intersectRayWithShape(sphere, ray, hitInfo);
        
        if (features.enableTextureMapping) {
            glm::vec3 bary = computeBarycentricCoord(v0Final.position, v1Final.position, v2Final.position, ray.origin + ray.t * ray.direction);
            glm::vec2 texCoord = interpolateTexCoord(v0Final.texCoord, v1Final.texCoord, v2Final.texCoord, bary);
            hitInfo.texCoord = texCoord;
            if (hitInfo.material.kdTexture) {
                hitInfo.material.kd = acquireTexel(*hitInfo.material.kdTexture.get(), texCoord, features);
            }
        }

        if (features.enableNormalInterp == true) {
            glm::vec3 bary = computeBarycentricCoord(v0Final.position, v1Final.position, v2Final.position, ray.origin + ray.t * ray.direction);
            glm::vec3 interpolatedNormal = interpolateNormal(v0Final.normal, v1Final.normal, v2Final.normal, bary);
            Ray n0, n1, n2, nIterpolated;
            n0.origin = v0Final.position;
            n0.direction = v0Final.normal;
            n0.t = 1.0f;
            n1.origin = v1Final.position;
            n1.direction = v1Final.normal;
            n1.t = 1.0f;
            n2.origin = v2Final.position;
            n2.direction = v2Final.normal;
            n2.t = 1.0f;
            nIterpolated.origin = ray.origin + ray.t * ray.direction;
            nIterpolated.direction = interpolatedNormal;
            if (glm::dot(glm::normalize(ray.direction), interpolatedNormal) > 0)
                nIterpolated.direction = -1.0f * interpolatedNormal;
            nIterpolated.t = 1.0f;
            drawRay(n0);
            drawRay(n1);
            drawRay(n2);
            drawRay(nIterpolated);
        }
        //Finish timing
        #ifdef traversal_debug_time
        auto t2 = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double, std::milli> delta = t2 - t1;
        std::cout << "Time elapsed: " << delta.count() << " ms\n";
        #endif
        return hit;
        

    }

    
}
