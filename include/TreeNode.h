#include <iostream>
#include <vector>

using namespace std;

class TreeNode {
public:
    size_t id;
    TreeNode* parent;
    vector<TreeNode*> children;

    TreeNode(size_t id) : id(id), parent(nullptr) {}

    void addChild(TreeNode* child);
    void removeChild(TreeNode* child);
    TreeNode* findChildById(size_t id);
    vector<TreeNode*> getChildren();

    vector<TreeNode*> getLeaves();

private:
    void getLeavesHelper(const TreeNode* node, vector<TreeNode*>& leaves) const;
};

