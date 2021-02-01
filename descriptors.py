# descriptors.py
class Verbose_attribute():
    def __get__(self, obj, type=None) -> object:
        print("accessing the attribute to get the value")
        return 42
    def __set__(self, obj, value) -> None:
        print("accessing the attribute to set the value")
        #this is a read-only descriptor
        raise AttributeError("Cannot change the value")

class Foo():
    attribute1 = Verbose_attribute()
    #print("meh", attribute1.__dict__)

my_foo_object = Foo()
#print(my_foo_object)
x = my_foo_object.attribute1
print(x)