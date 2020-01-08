class BidirectionalNode:
    def __init__(self, value=None, pred=None, succ=None):
        self.pred = pred
        self.succ = succ
        self.value = value

    def __repr__(self):
        return repr(self.value)


class BidirectionalList:
    def __init__(self, items=None):
        self.head = None
        self.tail = None

        if items:
            for item in items:
                self.append(item)

    def prepend(self, value):
        # TODO Prepend multiple if value is list
        new_head = BidirectionalNode(value=value, succ=self.head)
        if self.head:
            self.head.pred = new_head
        if not self.tail:
            self.tail = new_head
        self.head = new_head
        return new_head

    def append(self, value):
        # TODO Append multiple if value is list
        new_tail = BidirectionalNode(value=value, pred=self.tail)
        if self.tail:
            self.tail.succ = new_tail
        if not self.head:
            self.head = new_tail
        self.tail = new_tail
        return new_tail

    def insert_after(self, after, value):
        node = BidirectionalNode(value=value, pred=after, succ=after.succ)
        if after.succ:
            after.succ.pred = node
        after.succ = node
        if self.tail is after:
            self.tail = node
        return node

    def insert_before(self, before, value):
        node = BidirectionalNode(value=value, pred=before.pred, succ=before)
        if before.pred:
            before.pred.succ = node
        before.pred = node
        if self.head is before:
            self.head = node
        return node

    def remove(self, node):
        if node.pred:
            node.pred.succ = node.succ
        if node.succ:
            node.succ.pred = node.pred
        if node is self.head:
            self.head = node.succ
        if node is self.tail:
            self.tail = node.pred
        node.pred = None
        node.succ = None

    def __getitem__(self, idx):
        if idx >= 0:
            curr = self.head
            while curr:
                if idx == 0:
                    return curr
                curr = curr.succ
                idx -= 1
        else:
            idx += 1
            curr = self.tail
            while curr:
                if idx == 0:
                    return curr
                curr = curr.pred
                idx += 1
        raise KeyError

    def __len__(self):
        count = 0
        curr = self.head
        while curr:
            count += 1
            curr = curr.succ
        return count

    def __iter__(self):
        curr = self.head
        while curr:
            yield curr
            curr = curr.succ

    def __repr__(self):
        return repr(list(self))
