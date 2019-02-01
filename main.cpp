#include <algorithm>
#include <iostream>
#include <stdio.h>
#include <random>
#include <math.h>
#include <time.h>
#include <vector>
#include <queue>
#include <set>

using namespace std;

#define MAX 16
#define INF 1000000000
#define EPS 1e-9
#define MASK 1000000
#define WALL -10000000
#define OP_NUM 100000

template<typename T> inline void Maxi(T &x, T y) {if(x < y) x = y; }
template<typename T> inline void Mini(T &x, T y) {if(x > y) x = y; }

struct Point
{
    int x, y;
    Point operator + (Point p) { return {x + p.x, y + p.y}; }
    Point operator - (Point p) { return {x - p.x, y - p.y}; }
    Point operator * (int k) { return {x * k, y * k}; }
    bool operator == (Point p) { return x == p.x && y == p.y; }
};

int distance(Point p)
{
    return abs(p.x) + abs(p.y);
}

struct Cell
{
    int x, y, val, dis;
    bool operator < (Cell c)
    {
    	if(dis <= 1 && c.dis > 1) return 1;
    	if(dis > 1 && c.dis <= 1) return 0;
    	return val < c.val;
    }
};

Point dir[4] = {{-1, 0},
                { 0, 1},
                { 1, 0},
                { 0,-1}};

class Bot
{
public:
    Bot() {}
    ~Bot() {}

    Point curPos;
    Point nextPos;
    int nTrash;
    bool isMask;
    bool isCrash;

    int maxVal[MAX][MAX];
    int dist[MAX][MAX];
    int dirRoot[MAX][MAX];

    Bot clone()
    {
        Bot bot;
        bot.curPos = {curPos.x,curPos.y};
        bot.nextPos = {nextPos.x,nextPos.y};
        bot.nTrash = nTrash;
        bot.isMask = isMask;
        bot.isCrash = isCrash;
        return bot;
    }

    void update(int board[MAX][MAX], int nRow, int nCol)
    {
        for(int i = 1; i <= nRow; ++i) {
            for(int j = 1; j <= nCol; ++j) {
                dist[i][j] = -1;
                maxVal[i][j] = -INF;
            }
        }
        queue<Point> Qu;
        Qu.push(curPos);
        dist[curPos.x][curPos.y] = 0;
        maxVal[curPos.x][curPos.y] = 0;
        while(!Qu.empty())
        {
            Point cur = Qu.front();
            Qu.pop();
            for(int i = 0; i < 4; ++i)
            {
                Point next = cur + dir[i];
                if(next.x < 1 || next.x > nRow || next.y < 1 || next.y > nCol) continue;
                if(!isMask && board[next.x][next.y] == WALL) continue;
                int val = maxVal[cur.x][cur.y] + (isMask ? max(0,board[next.x][next.y]) : board[next.x][next.y]);
                if(maxVal[next.x][next.y] < val)
                {
                    maxVal[next.x][next.y] = val;
                    if(cur == curPos) dirRoot[next.x][next.y] = i;
                    else dirRoot[next.x][next.y] = dirRoot[cur.x][cur.y];
                }
                if(dist[next.x][next.y] >= 0) continue;
                dist[next.x][next.y] = dist[cur.x][cur.y] + 1;
                Qu.push(next);
            }
        }
    }

    void present()
    {
        cerr<<curPos.x<<' '<<curPos.y<<'\n';
        cerr<<nTrash<<' '<<isMask<<' '<<isCrash<<'\n';
    }
};

void fread(FILE *f, int &x)
{
    char ch = 0;
    while(!(('0' <= ch && ch <= '9') || ch == 'M' || ch == 'W'))
    {
        if(fscanf(f,"%c",&ch) == -1) break;
    }
    if(ch == 'M') x = MASK;
    else if(ch == 'W') x = WALL;
    else {
        x = 0;
        do {
            x = 10 * x + ch - '0';
            if(fscanf(f,"%c",&ch) == -1) break;
        } while('0' <= ch && ch <= '9');
    }
}

class GameState
{
public:
    GameState() {}
    ~GameState() {}

    void init()
    {
        FILE *f = fopen("MAP.inp", "r");
        if(f == NULL) return;
        fscanf(f,"%d%d%d", &nRow, &nCol, &nMove);
        for(int i = 0; i < 2; ++i) {
            fscanf(f, "%d%d", &bot[i].curPos.x, &bot[i].curPos.y);
        }
        fscanf(f, "%d%d", &bot[0].nTrash, (int*)&bot[0].isMask);
        for(int i = 1; i <= nRow; ++i) {
            for(int j = 1; j <= nCol; ++j) {
                // TODO: read string
                fread(f,board[i][j]);
            }
        }
        fclose(f);
    }

    GameState clone()
    {
        GameState st;
        st.nMove = nMove;
        st.nRow = nRow;
        st.nCol = nCol;
        for(int i = 1; i <= nRow; ++i) {
            for(int j = 1; j <= nCol; ++j) st.board[i][j] = board[i][j];
        }
        st.bot[0] = bot[0].clone();
        st.bot[1] = bot[1].clone();
        return st;
    }

    void doMove(int move[])
    {
        bot[0].nextPos = bot[0].curPos + dir[move[0]];
        bot[1].nextPos = bot[1].curPos + dir[move[1]];
        nMove--;
        
        for(int i = 0; i < 2; ++i)
        {
            int x = bot[i].nextPos.x;
            int y = bot[i].nextPos.y;
            if(board[x][y] == MASK) bot[i].isMask = 1;
            else bot[i].nTrash += bot[i].isMask ? max(0,board[x][y]) : board[x][y];
            if(board[x][y] > 0) board[x][y] = 0;
            bot[i].curPos = bot[i].nextPos;
        }

        if(bot[0].nextPos == bot[1].nextPos)
        {
            bot[0].isCrash = bot[1].isCrash = 1;
            return;
        }
    }

    multiset<int> getMove(int player)
    {
        multiset<int> move;
        if(nMove == 0 || bot[player].isCrash) return move;

        vector<Cell> St;
        St.clear();
        bot[player].update(board, nRow, nCol);
        for(int i = 1; i <= nRow; ++i) {
            for(int j = 1; j <= nCol; ++j) {
                int plr = distance(bot[player].curPos - (Point){i,j});
                int opn = distance(bot[player^1].curPos - (Point){i,j});
                if(opn == plr && bot[player].nTrash < bot[player^1].nTrash) continue;
                if(opn-plr < 0) continue;
                St.push_back({i,j,-bot[player].maxVal[i][j],plr});
            }
        }
        sort(St.begin(),St.end());
        if(St[0].val == 0)
        {
        	for(int i = 1; i <= nRow; ++i) {
	            for(int j = 1; j <= nCol; ++j) {
	                St.push_back({i,j,-bot[player].maxVal[i][j],0});
	            }
	        }
	        sort(St.begin(),St.end());
        }
        for(int i = 0; i < min((int)St.size(),32); ++i) {
            int x = St[i].x;
            int y = St[i].y;
            move.insert(bot[player].dirRoot[x][y]);
        }
        St.clear();
        return move;
    }

    int getResult(int player)
    {
        return bot[player].nTrash;
    }

    void present()
    {
        cerr<<"Bot 0:\n";
        bot[0].present();
        cerr<<"Bot 1:\n";
        bot[1].present();
        cerr<<nMove<<'\n';
        for(int i = 1; i <= nRow; ++i) {
            for(int j = 1; j <= nCol; ++j)
            {
                if(bot[0].curPos.x == i && bot[0].curPos.y == j) 
                    cerr<<"\033[1;4;31mA\033[0m ";
                else if(bot[1].curPos.x == i && bot[1].curPos.y == j) 
                    cerr<<"\033[1;4;34mB\033[0m ";
                else
                {
                    if(board[i][j] == MASK) cerr<<"M ";
                    else if(board[i][j] == WALL) cerr<<"W ";
                    else cerr<<board[i][j]<<' ';
                }
            }
            cerr<<'\n';
        }
    }

    int nMove;
    int nRow, nCol;
    int board[MAX][MAX];
    Bot bot[2];
};


class Node
{
public:
    Node(GameState _state, int _move[] = NULL, Node *_parentNode = NULL)
    {
        parentNode = _parentNode;
        childNodes.clear();
        if(_move != NULL)
        {
            move[0] = _move[0];
            move[1] = _move[1];
        }
        
        untriedMove[0] = _state.getMove(0);
        untriedMove[1] = _state.getMove(1);
    }
    ~Node() {}

    double UCB(Node *node, int player)
    {
        return node->scores[player] / node->visits + sqrt(2.0 * log(visits) / node->visits);
    }

    Node* UCTSelectChild()
    {
        Node *node = childNodes[0];
        for(int i = 1; i < (int)childNodes.size(); ++i) {
            if( UCB(node,0) + EPS < UCB(childNodes[i],0)
            || (fabs(UCB(node,0) - UCB(childNodes[i],0) < EPS)
                && UCB(node,1) + EPS < UCB(childNodes[i],1)) ) node = childNodes[i];
        }
        return node;
    }

    Node* addChild(int move[], GameState s)
    {
        Node *node = new Node(s, move, this);
        untriedMove[0].erase(untriedMove[0].find(move[0]));
        untriedMove[1].erase(untriedMove[1].find(move[1]));
        childNodes.push_back(node);
        return node;
    }

    void update(int r0, int r1)
    {
        scores[0] += r0;
        scores[1] += r1;
        visits += 1;
    }

    int move[2];
    Node *parentNode;
    vector<Node*> childNodes;
    double scores[2];
    double visits;
    multiset<int> untriedMove[2];
};

int random(multiset<int> s)
{
    int r = rand() % (int)s.size();
    multiset<int>::iterator it = s.begin();
    while(r)
    {
        it++;
        r--;
    }
    return *it;
}

double result[2];

int* UCT(GameState rootState, int itermax)
{
    Node *rootNode = new Node(rootState);
    for(int i = 0; i < itermax; ++i)
    {
        Node *node = rootNode;
        GameState state = rootState.clone();

        while((node->untriedMove[0].empty() || node->untriedMove[1].empty()) && !node->childNodes.empty())
        {
            node = node->UCTSelectChild();
            state.doMove(node->move);
        }

        if(!node->untriedMove[0].empty() && !node->untriedMove[1].empty())
        {
            int move[2];
            move[0] = random(node->untriedMove[0]);
            move[1] = random(node->untriedMove[1]);
            state.doMove(move);
            node = node->addChild(move,state);
        }

        while(!state.getMove(0).empty() && !state.getMove(1).empty())
        {
            int move[2];
            move[0] = random(state.getMove(0));
            move[1] = random(state.getMove(1));
            state.doMove(move);
        }

        while(node != NULL)
        {
            node->update(state.getResult(0),state.getResult(1));
            node = node->parentNode;
        }
    }
    
    Node *node = rootNode->childNodes[0];
    for(int i = 1; i < (int)rootNode->childNodes.size(); ++i) {
        if(node->visits < rootNode->childNodes[i]->visits) {
            node = rootNode->childNodes[i];
        }
    }
    result[0] = node->UCB(node,0);
    result[1] = node->UCB(node,1);
    return node->move;
}

void playMove(Point pos)
{
    FILE *f = fopen("MOVE.out", "w");
    fprintf(f, "%d %d\n", pos.x, pos.y);
    fclose(f);
}

bool vst[MAX][MAX];

void firstPlay(GameState *state)
{
    vector<Point> Mask;
    for(int x = 1; x <= state->nRow; ++x) {
        for(int y = 1; y <= state->nCol; ++y) {
            if(state->board[x][y] == MASK) {
            	Mask.push_back((Point){x,y});
            }
        }
    }
	GameState st = state->clone();
	st.bot[0].curPos = Mask[0];
	st.bot[1].curPos = Mask[1];
	st.nMove -= 2;
	UCT(st, OP_NUM / (st.nMove * st.nRow * st.nCol));
	int pl;
	if(result[0] > result[1]) pl = 0;
	else pl = 1;
	queue<Point> Qu;
	Qu.push(Mask[pl]);
	vst[Mask[pl].x][Mask[pl].y] = 1;
	while(!Qu.empty())
	{
		Point pos = Qu.front();
		Qu.pop();
		for(int i = 0; i < 4; ++i)
		{
			Point next = pos + dir[i];
			if(vst[next.x][next.y]) continue;
            if(state->board[next.x][next.y] == WALL) continue;
			if(state->board[next.x][next.y] == 0)
			{
				playMove(next);
				return;	
			}
			Qu.push(next);
			vst[next.x][next.y] = 1;
		}
	}	
}

void UCTPlayGame()
{
    srand(time(NULL));
    GameState *state = new GameState();
    state->init();

    if(state->bot[0].curPos == (Point){0,0}) firstPlay(state);
    else
    {
        int* move = UCT(*state, OP_NUM / (state->nMove * state->nRow * state->nCol));
        state->doMove(move);
        playMove(state->bot[0].curPos);    
    }
}

int main()
{   
    UCTPlayGame();
    return 0;
}