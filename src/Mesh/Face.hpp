#ifndef FACE_HPP
#define FACE_HPP

#include <vector>
#include <memory>

#include "../Vars.hpp"

class Face
{
    public:
        enum faceType{GENERAL, TRIANGULAR, QUADRILATERAL};

        Face() : type(GENERAL) {};
        Face(std::vector<int> nodesIdx) : nodesIndex(nodesIdx), type(GENERAL) {};
        Face(std::vector<int> nodesIdx, faceType fType) : nodesIndex(nodesIdx), type(fType) {};

        void update(const std::vector<Vars<3>>& nodeList);
        bool check() const;
        bool equal(const Face& compFace) const;

        void reverseOrientation();

        virtual ~Face();

        std::vector<int> nodesIndex;

        double area;
        Vars<3> normalVector;
        Vars<3> midpoint;

    protected:
        const int type;

        Vars<3> calculateNormalVector(const std::vector<Vars<3>>& nodeList);
        Vars<3> calculateMidpoint(const std::vector<Vars<3>>& nodeList) const;

};

#endif // FACE_HPP