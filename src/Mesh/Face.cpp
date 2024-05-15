#include <cmath>
#include <algorithm>
#include "Face.hpp"


void Face::update(const std::vector<Vars<3>>& nodeList)
{
    midpoint = calculateMidpoint(nodeList);

    Vars<3> normalVectorScale = calculateNormalVector(nodeList);

    area = norm2(normalVectorScale);
    normalVector = normalVectorScale/area;
}

bool Face::check() const
{
    if(std::adjacent_find(nodesIndex.begin(), nodesIndex.end()) != nodesIndex.end())
    {
        return true;
    }

    return false;
}

bool Face::equal(const Face& compFace) const
{
    bool isEqual = true;

    for (auto & compNodeIndex : compFace.nodesIndex)
    {
        if(std::find(nodesIndex.begin(), nodesIndex.end(), compNodeIndex) == nodesIndex.end())
        {
            isEqual = false;
        }
    }
    
    return isEqual;
}

void Face::reverseOrientation()
{
    std::reverse(nodesIndex.begin(), nodesIndex.end());
    //update();
}

Vars<3> Face::calculateNormalVector(const std::vector<Vars<3>>& nodeList)
{
    Vars<3> surface = Vars<3>();

    int i;
    for (i = 0; i < nodesIndex.size() - 1; i++)
    {
        Vars<3> auxSurface = cross(nodeList[nodesIndex[i+1]] - nodeList[nodesIndex[i]], (midpoint - nodeList[nodesIndex[i]])/2.0);
        surface = surface + auxSurface;
    }
    
    Vars<3> auxSurface = cross(nodeList[nodesIndex[0]] - nodeList[nodesIndex[i]], (midpoint - nodeList[nodesIndex[i]])/2.0);
    surface = surface + auxSurface;

    return surface;
}

Vars<3> Face::calculateMidpoint(const std::vector<Vars<3>>& nodeList) const
{
    Vars<3> aux = Vars<3>();

    int i;
    for (i = 0; i < nodesIndex.size(); i++)
    {
        aux = aux + nodeList[nodesIndex[i]];
    }
    
    return aux / ((double) i);
}

Face::~Face()
{
    
}