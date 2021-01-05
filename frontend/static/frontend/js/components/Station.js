import React, { Component, Suspense } from "react";
import "../App.css";
import Loading from "./LoadingScreen";
import Routes from "./Routes";
import { BrowserRouter as Router } from "react-router-dom";
import { Row, Col, Card } from "@geist-ui/react";

export default class Station extends Component {
  constructor(props) {
    super(props);
    this.state = {
      res: [],
      mode: "",
    };
  }

  responseCatcher = (res, mode) => {
    this.setState({ res, mode });
  };

  render() {
    return (
      <Row justify="center">
        <Col span={18}>
          <Card>
            <Suspense fallback={<Loading />}>
              <Router>
                <Routes responseCatcher={this.responseCatcher} {...this.state} />
              </Router>
            </Suspense>
          </Card>
        </Col>
      </Row>
    );
  }
}
