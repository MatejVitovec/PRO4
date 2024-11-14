#ifndef VOLFIELD_H
#define VOLFIELD_H

#include "Field.hpp"
#include "Mesh/Mesh.hpp"

template <typename T>
class VolField : public Field<T>
{
    public:
        using Field<T>::operator+=;
        using Field<T>::operator-=;

        VolField() = delete;
        VolField(const Field<T>& field) : Field<T>(field), boundaryData() {}
        VolField(Mesh mesh_) : Field<T>(mesh_.getCellsSize()), boundaryData(mesh_.getBoundarySize()) {}
        VolField(Mesh mesh_, T def) : Field<T>(mesh_.getCellsSize(), def), boundaryData(mesh_.getBoundarySize()) {}

        virtual ~VolField() {}

        void initBoundary(int i)
        {
            boundaryData = std::vector<std::vector<T>>(i);
        }

        const std::vector<T>& boundary(int i) const
        {
            if (i >= boundaryData.size())
            {
                std::cout << "chyba pristupu do neuexistijiciho pole boundary values"
            }
            return boundaryData[i];
        }

        std::vector<T>& boundary(int i)
        {
            if (i >= boundaryData.size())
            {
                boundaryData.resize(index + 1);
            }
            return boundaryData[i];
        }

    private:
        std::vector<std::vector<T>> boundaryData;
};

#endif // VOLFIELD_H