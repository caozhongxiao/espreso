#ifndef ISAVABLEOBJECT_H
#define ISAVABLEOBJECT_H

namespace espreso {

    class ISavableObject
    {
    public:
        virtual ~ISavableObject() {};

        virtual void save() = 0;
    };

}

#endif // ISAVABLEOBJECT_H
