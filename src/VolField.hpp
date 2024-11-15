#ifndef VOLFIELD_H
#define VOLFIELD_H

#include <iostream>

#include "Field.hpp"
#include "Mesh/Mesh.hpp"

template <typename T>
class VolField : public Field<T>
{
    public:
        using Field<T>::operator+=;
        using Field<T>::operator-=;

        VolField() = delete;
        //VolField(const Field<T>& field) : Field<T>(field), boundaryData() {}
        VolField(Mesh mesh_) : Field<T>(mesh_.getCellsSize()), boundaryData(mesh_.getBoundarySize()) {}
        VolField(Mesh mesh_, T def) : Field<T>(mesh_.getCellsSize(), def), boundaryData(mesh_.getBoundarySize()) {}

        virtual ~VolField() {}

        void initBoundary(int i)
        {
            boundaryData = std::vector<std::vector<T>>(i);
        }

        size_t boundarySize() const
        {
            return boundaryData.size();
        }

        const std::vector<T>& boundary(int i) const
        {
            if (i >= boundaryData.size())
            {
                std::cout << "chyba pristupu do neuexistijiciho pole boundary values" << std::endl;
            }
            return boundaryData[i];
        }

        std::vector<T>& boundary(int i)
        {
            if (i >= boundaryData.size())
            {
                boundaryData.resize(i + 1);
            }
            return boundaryData[i];
        }

        const std::vector<std::vector<T>>& getBoundaryData() const
        {
            return boundaryData;
        }

        std::vector<std::vector<T>>& getBoundaryData()
        {
            return boundaryData;
        }

        void update(const Field<T>& w)
        {
            Field<T>::data = w.getData(); //predelat na friend
        }

    private:
        std::vector<std::vector<T>> boundaryData;
};

#endif // VOLFIELD_H