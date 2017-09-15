#ifndef IVALIDATABLEOBJECT_H
#define IVALIDATABLEOBJECT_H

namespace espreso
{

    class IValidatableObject
    {
    public:
        virtual ~IValidatableObject() {}

        virtual bool isValid() = 0;
        virtual QString errorMessage() = 0;
    };

}

#endif // IVALIDATABLEOBJECT_H
