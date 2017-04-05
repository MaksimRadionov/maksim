#include <string>
#include <stdio.h>
#include <cstdlib>
#include <time.h>
#include <set>
#include <iostream>
using namespace std;

struct AVLNode {
    AVLNode *left, *right; 
    string key;
	int value; //добавили
    int height;
	AVLNode(string const &_key, int _value) : key(_key), value(_value)
	{
       left = right = NULL;
       height = 1;
    }
    int getBalanceFactor() const {
        int r = right == NULL ? 0 : right->height;
        int l = left == NULL ? 0 : left->height;
        return r - l;
    }
    void fix()//определяет высоту узла по высотам детей 
	{
        int r = right == NULL ? 0 : right->height;
        int l = left == NULL ? 0 : left->height;
        height = (r > l ? r : l) + 1;
    }

     AVLNode *insert(string const &_key, int _value) 
	 {
         if (_key < key)
            left = left == NULL ? new AVLNode(_key,_value) : left->insert(_key,_value);
         else
            right = right == NULL ? new AVLNode(_key, _value) : right->insert(_key,_value);
         return balance();
     }
     AVLNode *findMinimum() {
        return left != NULL ? left->findMinimum() : this;
     }   

     AVLNode *removeMinimum() {
        if (left == NULL) return right;
        left = left->removeMinimum();
        return balance();
     }

     static AVLNode *remove(AVLNode *p, string const &_key) {
        if (p == NULL) return NULL;
        if (_key < p->key) {
           p->left = remove(p->left, _key);
           return p->balance();
        } else if (_key > p->key) {
           p->right = remove(p->right, _key);
           return p->balance();
        } else {
           AVLNode *l = p->left;
           AVLNode *r = p->right;
           delete p;
           if (r == NULL) return l;
           AVLNode *min = r->findMinimum();
           min->right = r->removeMinimum();
           min->left = l;
           return min->balance();
        }
    }
	 AVLNode* balance()
	 {
		 fix();
		 switch (getBalanceFactor()){
		 case -2:
			 if (left->getBalanceFactor() > 0)
				 left = left->rotateLeft();
			 return rotateRight();
		 case 2:
			 if (right->getBalanceFactor() < 0)
				 right = right->rotateRight();
			 return rotateLeft();
		 default:return this;
		 }
	 }
	 AVLNode * rotateLeft(){
		 AVLNode * t = right;
		 right = t->left;
		 t->left = this;
		 fix();
		 t->fix();
		 return t;

	 }
	 AVLNode * rotateRight(){
		 AVLNode * t = left;
		 left = t->right;
		 t->right = this;
		 fix();
		 t->fix();
		 return t;

	 }
	 void print(int indent)
	 {
		 if (right != NULL)
			 right->print(indent+1);
		 for (int i = 0; i < indent; i++)
		 {
			 printf("    ");
		 }
		printf("%s\n", key.c_str());
		if (left != NULL)
		left->print(indent + 1);
		
		 
	 }
	 AVLNode* find(string const &_key)
	 {
		// if (this == NULL) return NULL; то есть можно как бы обратиться NULL->find()
		 if (key == _key)return this;
		// if (key < _key)return left->find(_key);
		 if (key < _key)return left == NULL ? NULL : left->find(_key);
		 return right == NULL ? NULL : right->find(_key);

	 }

};
//моджифицировать вставку, если ключ уже есть то не втсавялть
struct AVLTree {
    AVLNode *root;
    AVLTree() {
        root = NULL;
    }
	~AVLTree()
	{
		remove(root);
	}

    void print() const 
	{
        if (root != NULL) root->print(0);
    }

    bool insert(string const &_key) {
        if (root == NULL) root = new AVLNode(_key);
        else root = root->insert(_key);
        return true;
    }

    bool remove(string const &_key) {
        root = AVLNode::remove(root, _key);
        return true;
    }
	bool find(string const &_key, int val)
	{
		if (root == NULL)return false;
		AVLNode*f = root->find(_key);
		if (f == NULL)return false;//return f==NUll?false:val=f->value,true; (оператор запятая, 
		val = f->value;
		return true;
	}
private:
	void remove(AVLNode * t)
	{
		if (t == NULL) return;
		remove(t->left);
		remove(t->right);
		delete t;
	}
	/*можно пергрузить оператор [] :
	int & AVLTree::operator[](string const &_key){      & чтобы можно было делать так: t["Париж"]=15000
	root->insert(_key,0);
	return find(_key).value
	*/
};
    
int main() 
{
    AVLTree t;
  /*  t.insert("abra");
    t.insert("cadabra");
    t.insert("foo");
    t.insert("bar");
    t.print();
    t.remove("cadabra");
    t.print();*/
	t.insert("1");
	t.insert("2");
	t.insert("3");
	t.insert("4");
	t.print();
	cout << endl;
	t.remove("1");
	t.print();
	cout << endl;
	int n;
	cin >> n;
}
