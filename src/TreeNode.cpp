//
// Created by Henry Hollis on 3/20/24.
//
#include <TreeNode.h>
void TreeNode::addChild(TreeNode* child) {
    children.push_back(child);
    child->parent = this;
}
void TreeNode::removeChild(TreeNode* child) {
    children.erase(std::remove(children.begin(), children.end(), child), children.end());
    child->parent = nullptr;
}

TreeNode* TreeNode::findChildById(size_t id) {
    for (auto* child : children) {
        if (child->id == id) {
            return child;
        }
    }
    return nullptr;
}

vector<TreeNode*> TreeNode::getChildren() {
        return children;
}

vector<TreeNode*> TreeNode::getLeaves() {
    vector<TreeNode*> leaves;
    getLeavesHelper(this, leaves);
    return leaves;
}

void TreeNode::getLeavesHelper(const TreeNode* node, vector<TreeNode*>& leaves)const {
    if (node->children.empty()) {
        leaves.push_back(const_cast<TreeNode*>(node));
    } else {
        for (const auto& child : node->children) {
            getLeavesHelper(child, leaves);
        }
    }
}

TreeNode* searchLeaves(const vector<TreeNode*>& leaves, size_t id) {
    for (TreeNode* leaf : leaves) {
        if (leaf->id == id) {
            return leaf;
        }
    }
    return nullptr;
}

vector<TreeNode*> move_node(vector<TreeNode*>& leaves, size_t from_node_id, size_t  to_node_id, size_t childID) {
    TreeNode* from_node = searchLeaves(leaves, from_node_id);
    TreeNode* to_node = searchLeaves(leaves, to_node_id);
    if (from_node && to_node){
        TreeNode* childToMove = from_node->findChildById(childID);
        if (childToMove) {
            // Remove the child branch from its current parent
            childToMove->parent->removeChild(childToMove);

            // Add the child branch as a child of toNode
            to_node->addChild(childToMove);

            // Check if from_node is now childless and remove it from leaves if needed
            if (from_node->children.size() == 0) {
                auto it = std::find(leaves.begin(), leaves.end(), from_node);
                if (it != leaves.end()) {
                    leaves.erase(it);
                }
            }

        }else
            cerr<<"Child with ID: "<< childID << " not found under node with ID: "<< from_node_id << "."<<endl;

    }else
        cerr<<"One or both nodeIDs provided were not found in vector provided."<<endl;
    return leaves;

}

void printTree(const vector<TreeNode*>& leaves, int depth = 0) {
    for (const auto& node : leaves) {
        for (int i = 0; i < depth; ++i) {
            cout << "  ";
        }
        if (node->parent) {
            cout << "|-- " << node->id << " (Parent: " << node->parent->id << ")" << endl;
        } else {
            cout << "|-- " << node->id << " (Root)" << endl;
        }
        printTree(node->getChildren(), depth + 1);
    }
}

vector<TreeNode*> mergeNodes(vector<TreeNode*>& leaves, size_t id1, size_t  id2, size_t parentID) {
    TreeNode* node1 = searchLeaves(leaves, id1);
    TreeNode* node2 = searchLeaves(leaves, id2);
    if (node1 && node2) {
        TreeNode *parent = new TreeNode(parentID);  // Placeholder ID for internal nodes
        parent->addChild(node1);
        if (node1 != node2)
            parent->addChild(node2);

        // Remove node1 and node2 from leaves vector
        leaves.erase(std::remove(leaves.begin(), leaves.end(), node1), leaves.end());
        if (node1 != node2)
            leaves.erase(std::remove(leaves.begin(), leaves.end(), node2), leaves.end());

        // Add the new parent node to the leaves vector
        leaves.push_back(parent);

    }else
        cerr<<"One or both nodeIDs provided were not found in vector provided."<<endl;

    return leaves;
}
