#include <iostream>

class Node{
  public:
    unsigned char direction;
    unsigned char permutation;
    Node* children[3];
    Node(); // Done
    Node(unsigned char, unsigned char); // Done
    ~Node(); // Done
    void addChild(unsigned char, unsigned char); // Done but unsecure
    void removeChildByDirection(unsigned char); // Done but unsecure
    int countChildren(); // Done
    void printTT(const std::string&, Node*, bool); // Done
    void cutBranch(std::string, Node*);
    void addBranch(std::string, std::string, Node*);
    unsigned char availablePermutation();
    std::string printForSerialization();
};

std::string serialize(Node*);
Node* beginUnserialize(std::string);
std::string generateAPlacement(std::string, Node*,int);
